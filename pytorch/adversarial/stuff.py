import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

import torchvision as tv
import torchvision.transforms as transforms
import torchvision.utils as utils

from IPython import display
import matplotlib.pyplot as plt


class LeNet5(nn.Module):
    def __init__(self):
        super(LeNet5, self).__init__()
        self.conv1 = nn.Conv2d(1,16,5)
        self.conv2 = nn.Conv2d(16,32,5)
        self.fc1 = nn.Linear(32*4*4,120)
        self.fc2 = nn.Linear(120,84)
        self.fc3 = nn.Linear(84,10)
        
    def forward(self, x):
        """
        Args:
          x of shape (batch_size, 1, 28, 28): Input images.
        
        Returns:
          y of shape (batch_size, 10): Outputs of the network.
        """
        bat = x.shape[0]
        x = F.max_pool2d(F.relu(self.conv1(x)),(2,2),stride=2)
        x = F.max_pool2d(F.relu(self.conv2(x)),(2,2),stride=2)
        x = x.view(bat,-1)#,32*4*4)
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        
        return self.fc3(x)

def compute_accuracy(net, testloader, device):
    net.eval()
    correct = 0
    total = 0
    with torch.no_grad():
        for images, labels in testloader:
            images, labels = images.to(device), labels.to(device)
            outputs = net(images)
            _, predicted = torch.max(outputs.data, 1)
            total += labels.size(0)
            correct += (predicted == labels).sum().item()
    return correct / total


def train_net(model: nn.Module, criteria, optimizer, epochs, tr_loader, ts_loader, device):
    loss_fn = criteria

    for i in range(epochs):
        model.train()
            
        for images,labels in tr_loader:
            images, labels = images.to(device), labels.to(device)
                
            yhat = model(images)
            loss = loss_fn(yhat,labels)
            loss.backward()
            optimizer.step()
            optimizer.zero_grad()
                
        print(compute_accuracy(model, ts_loader,device))
    

def plot_images(images, ncol=12, figsize=(8,8), cmap=plt.cm.Greys, clim=[0,1]):
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis('off')
    grid = utils.make_grid(images, nrow=ncol, padding=0, normalize=False).cpu()
    ax.imshow(grid[0], cmap=cmap, clim=clim)
    display.display(fig)
    plt.close(fig)

