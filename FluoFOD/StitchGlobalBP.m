function y=StitchGlobalBP(data)
s=size(data);
for i=1:s(2)
    vect=data(:,i);
    vect(vect==0)=[];
    c(i)=std(vect)/mean(vect);
    xi(i)=vect(1);
    yi(i)=vect(end);
end
low=median(c)-2*std(c);
high=median(c)+2*std(c);
% remove outliers in cov
rem=find((c<low)|(c>high));
display('Outliers removed prior to stitching')
numel(rem)
c(rem)=[];
xi(rem)=[];
yi(rem)=[];
data(:,rem)=[];
s=size(data);
%%
k=randperm(s(2));
for i=1:numel(k)
    list=[k(i)];
    xxi=xi;
    yyi=yi;
    current=k(i);
    % remove current start point
    xxi(k(i))=NaN;
    % while there are still unasigned pixels
    while sum(isnan(xxi))<numel(xxi)
        % find next fit
        dist=(yyi(current)-xxi).^2;
        idx=find(dist==min(dist),1,'first');
        list=[list,idx];
        % remove current endpoint
        yyi(current)=NaN;
        % update current
        current=idx;
        % remove current start point
        xxi(current)=NaN;
    end
    Rlist{i}=list;
    dataperm=data(:,list);
    dataperm(dataperm==0)=[];
    dataperm=dataperm(:);
    Rscore(i)=std(dataperm)/mean(dataperm); 
end
idx=find(Rscore==min(Rscore),1,'first');
data=data(:,Rlist{idx});
data(data==0)=[];
data=data(:);
y=data;    
