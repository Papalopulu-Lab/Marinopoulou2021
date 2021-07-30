function h=showfigure(x,m,raw,y1,md,M1,S1,S2,LLR,par2,i)
% adjust to fl intensity
m=m*S1+M1;
raw=raw*S1+M1;
y1=y1*S2*S1;
md=md*S2*S1;
%SHOWFIGURE Summary of this function goes here
%   Detailed explanation goes here
    % plots raw data with trend and detrended data
    h=figure();
    subplot(2,1,1)
    hold on
    plot(x,m,'r','LineWidth',2);
    plot(x,raw,'r');%,'+')
    xlabel('Time (hours)')
    ylabel('Fluorescence (a.u.)')
    str = sprintf('Raw fluorescent intensity cell %.0f',i);
    title(str,'fontweight','normal');
    subplot(2,1,2);    
    plot(x,y1);
    hold on
    plot(x,md,'b','LineWidth',2);
    xlabel('Time (hours)');
    ylabel('Detrended fluorescence (a.u.)');
    if LLR<0
    str = sprintf('LLR score %f',LLR);
    title(str);
    else
        newq=par2(2)/(2*pi*par2(1));
        if newq>10^8
            newq=Inf;
        end
    str = sprintf('LLR score %f, Period = %f, Quality = %f',LLR,(2*pi()/par2(2)),newq);
    title(str,'fontweight','normal');     
end

