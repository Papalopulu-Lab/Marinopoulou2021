function [LLRM,LL1,LL2,par1M,par2M,dataTOT]=getSynthLLRsimple(datamat,x,par1,par0,repeats)
dataTOT=[];
par1TOT=[];
% generate synthetic data
for i = 1:size(par1,1); % for each cell
    %     disp(i)
    y1=datamat(:,i);
    xcurr = x;
    xcurr(y1==0) = []; %deletes times from which no signal
    cov1=par1(i,2)*exp(-par1(i,1)*xcurr);
    Noise=par1(i,3);
    data=GetSynt(cov1,Noise,repeats,xcurr);
    datamp=zeros(repeats,numel(x));
    datamp(1:repeats,1:size(data,2))=data;
    dataTOT=[dataTOT;datamp];
    par1TOT=[par1TOT;repmat(par1(i,:),repeats,1)];
end
a0=mean(par1TOT(:,1));
parfor j=1:size(dataTOT,1)
    y1=dataTOT(j,:);
    xcurr=x; xcurr(y1==0)=[];
    y1(y1==0)=[];
    y1 = (y1 - mean(y1))/std(y1);
    [LLR,ll1,ll2, par1, par2] = fitsynthDataLLR(xcurr,y1',par1TOT(j,end),par1TOT(j,1));
    par1M(j,:)=par1;
    par2M(j,:) = par2;
    LLRM(j,:) = LLR;
    LL1(j,:)=ll1;
    LL2(j,:)=ll2;
end
end
%% functions needed
function data=GetSynt(cov1,Noise,repeats,x)
CovMatrix1 = zeros(length(x),length(x));
for j = 1:length(x)
    CovMatrix1(:,j) = circshift(cov1,[j-1,0]);
end
CovMatrix1 = CovMatrix1';
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';
MU = zeros(repeats,length(x));
Meas = diag((Noise^2).*ones(1,length(x)));
CVM1 = CVM1 + Meas;
SIGMA = CVM1; % change this to switch non-osc and osc
data= mvnrnd(MU,SIGMA);
end
