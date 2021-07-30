function [LLR,LL1,LL2,par1,par2,m] = fitsynthDataLLR(x,y1,noiseinit,varargin)
N=10; % number of runs
LLR = zeros(N,1);
par1vec = zeros(N,3);
par2vec = zeros(N,4);
avec=0.5+0.02*randn(100,1); avec(avec<0)=[]; avec=avec(1:N); 
bvec = 2 + 0.5*randn(N,1); 
if nargin>3
    a0=varargin{1};
    avec(6:end)=a0+a0/10*randn(5,1);
end
%% optimise OU model Kou=sigamsig*exp(-avec*t)+sigmanoise
for i = 1:N;
    samp = length(x);
    likfunc = @likGauss;
    par1mean = [log(avec(i)),log(var(y1))]; 
    % settings
    covfunc = @covOUa;
    hyp2.lik =log(noiseinit);
    hyp2.cov = par1mean;
    prior.lik ={{@priorDelta}};
    inf = {@infPrior,@infExact,prior};
    hyp2 = minimize(hyp2, @gp, -1000, inf, [], covfunc, likfunc, x, y1);
    % negative log likelihood- invert sign for log likelihood
    nlmlOU = gp(hyp2, @infExact, [], covfunc, likfunc, x, y1);
    % statistic used is 2log(ll)
    LLou(i)=-2*nlmlOU;%/numel(x);
    hypOU(i)=hyp2;
end
% get best OU model
[LL1,idx1] =max(LLou);
par1 = [exp(hypOU(idx1).cov), exp(hypOU(idx1).lik)];
%% fit the OU osc model
clear hyp2
for i=1:N
    covfunc = @covOUosca; 
    par2mean = [log(avec(i)),log(2*pi/bvec(i)),log(var(y1))];
    hyp2.lik =log(noiseinit);
    hyp2.cov = par2mean;
    prior.lik ={{@priorDelta}};
    inf = {@infPrior,@infExact,prior};
    hyp2 = minimize(hyp2, @gp, -1000, inf, [], covfunc, likfunc, x, y1); 
    nlmlOSC = gp(hyp2, @infExact, [], covfunc, likfunc, x, y1);
    LLosc(i)=-2*nlmlOSC;%/numel(x); 
    hypOUosc(i)=hyp2;
end
% choose best OU osc
[LL2,idx2] =max(LLosc);
par2 = [exp(hypOUosc(idx2).cov), exp(hypOUosc(idx2).lik)];
% generate fitting for the OUosc model onto data
xt=x;
[m,s2] = gp(hypOUosc(idx2),inf,[],@covOUosca,@likGauss,x,y1,xt);
%% Log-likelihood ratio 
if LL2>LL1
    LLR=(LL2-LL1)/numel(x)*100;
else
    LLR=0;
end
end
%%%%%%%%%%%%%%%%%

    