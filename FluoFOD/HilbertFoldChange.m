function y=HilbertFoldChange(dfit,dirname)
mkdir([dirname,'/Hilbert']);
for i=1:numel(dfit)
    x=dfit(i).x;
    m=dfit(i).m;
    M=dfit(i).M;
    fittedData=M+m;
    figure,subplot(2,1,1),plot(M,'r')
    hold on
    plot(fittedData,'b');
    title(['Cell ', num2str(i)]);
    % run hilbert on the slow trend
    tm=m-mean(m);
    hb=hilbert(tm);
    amp=abs(hb);
    ph=angle(hb);
    y(i).HilbertAmp=amp;
    y(i).HilbertPhase=ph;
    plot([M+mean(m)+amp],'b--');
    plot([M+mean(m)-amp],'b--');
    y(i).AbsoluteAmp=M+mean(m)+amp;
    subplot(2,1,2),
    plot(tm,'b');
    hold on
    plot(ph/max(ph)*max(tm),'g');
    legend('Detrended','Amplitude','Phase')
    % find minima points
    % find maxima points
    zcross=sign(ph(1:end-1))+sign(ph(2:end));
    tmp=find(zcross==0);
    list1=tmp(2:2:end);
    list2=tmp(1:2:end);
    % peak or trough
    if mean(tm(list1))>mean(tm(list2))
        hidx=list1;
        lidx=list2;
    else
        hidx=list2;
        lidx=list1;
    end
    % start from first peak- discard troughs before first peak
    lidx(lidx<hidx(1))=[];
    % pair the indices for high and low- start with high
    templidx=lidx;
    clear pair
    for j=1:numel(hidx)
        pair(j,1)=hidx(j);
        % index for the high levels
        if ~isempty(templidx)
            pair(j,2)=templidx(1);
            templidx(1)=[];
        else
            pair(j,2)=-1;
        end
    end
    % remove unpaired peaks and peaks too close together
    idx=find(pair(:,2)==-1);
    pair(idx,:)=[];
    diff=abs(pair(:,1)-pair(:,2));
    idx=find(diff<3); % was 5
    pair(idx,:)=[];
    plot(pair(:,1),tm(pair(:,1)),'r+');
    plot(pair(:,2),tm(pair(:,2)),'ko');
    % contains paired peaks (col1) with troughs (col2)
    y(i).PairedIdx=pair; 
    % calculate fold change from fitted signal
    fittedsig=m+M;
    fc1=fittedsig(pair(:,1))./fittedsig(pair(:,2));
    y(i).FoldChangeFit=fc1';
    % fold change from raw data
    subplot(2,1,1)
    for j=1:numel(fc1)
        text(pair(j,1),fittedsig(pair(j,1))+100,[num2str(round(fc1(j)*100)/100),'x']);
    end
    ylim=get(gca,'YLim');
    set(gca,'YLim',[ylim(1),ylim(2)+200]);
    subplot(2,1,1),legend('Trend','Model fit','Amplitude');
    subplot(2,1,2),legend('Detrended fit','Phase','Local peaks','Local troughs','Orientation','Horizontal','Location','NorthOutside');
    print(gcf(),[dirname,'/Hilbert/Cell ', num2str(i), 'Hilbert.png'],'-dpng');
end

%%
