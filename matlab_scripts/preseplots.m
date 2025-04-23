clear



nSamp=5000;
skew = 1.3;
kurt = skew^2 + 3.8;
pA=[0.05,0.1];
p0=[0.9,0.95];
rnge=8:.1:17;

m0 = 11;
mA = 15;

H0 = m0 + normrnd(0,.5,nSamp,1);
HA = mA + pearsrnd(0,1,-skew,kurt,nSamp,1);
q0=quantile(H0,p0);
qA=quantile(HA,pA);

map = brewermap(2,'Set1'); 
figure
histogram(H0,rnge,'FaceColor',map(2,:),'facealpha',.5,'edgecolor','none')
hold on
histogram(HA,rnge,'FaceColor',map(1,:),'facealpha',.5,'edgecolor','none')
box off
axis tight
%legend boxoff
xline(q0)
xline(qA)

HA = mA + normrnd(0,.5,nSamp,1);
H0 = m0 + pearsrnd(0,1,skew,kurt,nSamp,1);
q0=quantile(H0,p0);
qA=quantile(HA,pA);

map = brewermap(2,'Set1'); 
figure
histogram(H0,rnge,'FaceColor',map(2,:),'facealpha',.5,'edgecolor','none')
hold on
histogram(HA,rnge,'FaceColor',map(1,:),'facealpha',.5,'edgecolor','none')
box off
axis tight
%legend boxoff
xline(q0)
xline(qA)

% P(extreme inside 0.99 | samplesize=10) ~ 0.90
% P(extreme inside 0.95 | samplesize=10) ~ 0.60
