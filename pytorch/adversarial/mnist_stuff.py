from datetime import datetime

import torch
from torch import nn
from torch.nn import functional as F
from torch import optim
from torch.utils import data
import torchvision as tv


transforms: tv.transforms.Compose = tv.transforms.Compose([
        tv.transforms.ToTensor(),
        tv.transforms.Normalize([0.5], [0.5])])


class MnistModel(nn.Module):
    def __init__(self) -> None:
        super(MnistModel, self).__init__()

        self.block = nn.Sequential(
            nn.Conv2d(1, 32, 2),
            nn.MaxPool2d(2),
            nn.ReLU(),
            nn.Conv2d(32, 64, 2),
            nn.MaxPool2d(2),
            nn.ReLU(),
            nn.Conv2d(64, 128, 2),
            nn.ReLU(),
        )
        self.fc1 = nn.Linear(128 * 5**2 , 200)
        self.fc2 = nn.Linear(200, 10)
        self.relu = nn.ReLU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        out = self.block(x)
        out = out.view(-1,  128 * 5**2)
        out = self.relu(self.fc1(out))
        out = self.fc2(out)
        return F.log_softmax(out, 1)


def load_state(model: nn.Module, filename: str) -> nn.Module:
    state_dict = torch.load(filename)
    print("Loading state from: {}".format(filename))
    model.load_state_dict(state_dict)
    return model


def save_state(model: nn.Module, filename: str) -> None:
    torch.save(model.state_dict(), filename)
    print("Model saved to:", filename)


def now() -> str:
    return datetime.now().strftime('%b%d_%H-%M-%S')


def get_train_loader() -> data.DataLoader:
    train_set = tv.datasets.MNIST("/coursedata/", train=True, transform=transforms, download=True)
    return data.DataLoader(train_set, batch_size=32, shuffle=True)


def get_test_loader() -> data.DataLoader:
    test_set = tv.datasets.MNIST("/coursedata/", train=False, transform=transforms, download=True)
    return data.DataLoader(test_set, batch_size=32)


def train_mnist() -> nn.Module:
    train_loader: data.DataLoader = get_train_loader()
    test_loader: data.DataLoader = get_test_loader()

    model = MnistModel()

    optimizer = optim.Adam(model.parameters(), lr=0.005)
    criterion = F.cross_entropy

    print("Training model...")
    for epoch in range(10):
        print("Epoch {}/{}".format(epoch + 1, 10))

        running_loss = 0.0
        for i, (inputs, yreal) in enumerate(train_loader):

            optimizer.zero_grad()

            ypred = model(inputs)
            loss = criterion(ypred, yreal)
            loss.backward()
            optimizer.step()

            # book-keeping
            running_loss += loss.item()
            if i % 500 == (500 - 1):
                print('[%d/%d, %5d] loss: %.3f' % (epoch + 1, 10, i + 1, running_loss / 500))
                running_loss = 0.0

    test(model, test_loader, use_cuda=False)
    return model

def test(model: nn.Module, test_loader: data.DataLoader, use_cuda: bool) -> None:
    model.eval()

    if use_cuda:
        model = model.cuda()

    def test_average() -> None:
        correct = 0
        total = 0
        with torch.no_grad():
            for (inputs, yreal) in test_loader:
                if use_cuda:
                    inputs, yreal = inputs.cuda(), yreal.cuda()

                ypred = model(inputs)
                _, predicted = torch.max(ypred.data, 1)

                total += yreal.size(0)
                correct += (predicted == yreal).sum().item()

        accuracy = 100 * correct / total
        print("Accuracy of the network on the {} test images (average): {}".format(total, accuracy))

    def test_per_class() -> None:
        class_correct = list(0. for _ in range(10))
        class_total = list(0. for _ in range(10))
        total = 0

        with torch.no_grad():
            for (inputs, yreal) in test_loader:
                if use_cuda:
                    inputs, yreal = inputs.cuda(), yreal.cuda()

                total += yreal.size(0)

                ypred = model(inputs)
                _, predicted = torch.max(ypred, 1)
                c = (predicted == yreal).squeeze()
                for i in range(yreal.shape[0]):
                    label = yreal[i]
                    class_correct[label] += c[i].item()
                    class_total[label] += 1

        print("Accuracy of the network on the {} test images (per-class):".format(total))

        per_class_accuracy = {}
        for i in range(10):
            accuracy = 100 * class_correct[i] / (class_total[i] + 0.0001)
            per_class_accuracy[i] = accuracy
            print('Accuracy of %5s : %2d %%' % (
                i, accuracy))

    test_average()
    test_per_class()
