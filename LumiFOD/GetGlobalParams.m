function [posc,m,LL]=GetGlobalParams(data,time,len,varargin)
% len is detrending parameter
data=StitchGlobalBP(data);
dt=time(2)-time(1);
time=([1:numel(data)]-1)*dt;
time=time(:);
% reduce the length of global data to speed simulation up
% figure,subplot(2,1,1),plot(time,data)
% hold on
% disp('Performing global detrending..')
% try
%     data=data(1:1000);
%     time=time(1:1000);
% end
%% detrend the global data
if nargin>3
    m=varargin{1};
else
    [m,par0] = detrenddataNEW(data,time,len);
end
ddata=data-m;
%% MLE global estimate using OUosc covariance model
disp('Optimising global fit..')
N=1;
ell=2+0.01*randn(100,1); ell(ell<0)=[]; ell=ell(1:N);
y=ddata;
x=time;
pGlobal=zeros(N,4);
LL=zeros(N,1);
sf=var(y)
sn=0.10* var(y);
for i=1:N
    % cov model
    covfunc = @covOUosca;
    hyp.cov = log([ell(i); pi/2; sf]); % assumed period is 4h
    % noise term
    likfunc=@likGauss;
    hyp.lik =log(sn);
    % put prior on the noise
    prior.lik ={{@priorSmoothBox1,log(0.01*sf),log(0.15*sf),1000}};
    % inference
    inf = {@infPrior,@infExact,prior};
    % MLE
    hyp = minimize(hyp, @gp, -500, inf, [], covfunc, likfunc, x,y);
    nLL =gp(hyp, @infExact, [], covfunc, likfunc, x, y);
    LL(i)=-2*nLL;
    pGlobal(i,:)=exp([hyp.cov(1),hyp.cov(2),hyp.cov(3), hyp.lik]);
end
idx=find(LL==max(LL),1,'first');
posc=pGlobal(idx,:);
hyp.cov=log([posc(1),posc(2),posc(3)]);
hyp.lik=log(posc(end));
inf = {@infPrior,@infExact,[]};
[fm,sm] = gp(hyp,inf,[],@covOUosca,@likGauss,x,y,x);
ell0=posc(1)
sn0=posc(4)/posc(3)
SFG=posc(3)
p=2*pi/posc(2)