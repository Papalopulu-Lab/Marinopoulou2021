clear all, close all, clc
addpath(genpath('../GPML'))
startup
%% Step 1. Choose lengthscale to detrend (in hours) and other parameters
Lengthscale =7.5; % roughly 3x period expected (in hours) % 
DetrendParam = log(1./(2*Lengthscale.^2));
%% experiment details- this section can be duplicated to run more files, i.e. fnames{2}...
fnames{1}='20180915_Prolif.xls';
exptnames{1}='20180915_Prolif'; % custom name for the simulation file 
coldata{1}=[2:41];
bkgrdata{1}=[40:41];
strow{1}=1;% start row excludes rows with numbers in header
q = 0.05; % FDR recommended 1-5% 
for k=1:numel(fnames)
    dirname=exptnames{k};
    detrendData=[];
    modelFit=[];
    mkdir(dirname)
    num = xlsread(fnames{k});
    num(isnan(num)) = 0;
    data = num(strow{k}:end,coldata{k});
    bckgd=num(strow{k}:end,bkgrdata{k});
    % assigns time TIME ALWAYS COL 1
    time = num(strow{k}:end,1); % enter column of time
    time=(time*10)/60;%convert to hours
    %% remove cells that too short-optional
    crit=size(data,1)-sum(data==0);
    len=round(5.5/(time(2)-time(1))); % length cutoff is 5.5h
    rem=find(crit<len); % find the removed cells
    data(:,rem)=[];
    %set up matrices to fill
    par0TOT = zeros(size(data,2),3); 
    par1TOT = zeros(size(data,2),3);
    par2TOT = zeros(size(data,2),4);
    LLRM = zeros(size(data,2),1);
    PassList=zeros(size(data,2),1);
    LLRTOT=zeros(size(data,2),1);
    S1=[]; S2=[];
    dfit=struct('m',[],'M',[],'x',[]); % this stores the model fit
    %% Step 2. Estimate local variance of background-use as proxy for noise
    [stdev] = bckgdstdev(bckgd,time,-7);
    %% Step 3. Fit individual cells
    parfor i =1:size(data,2)
        i
        y1 =data(:,i);
        x = time;
        x(y1==0) = []; %deletes times from which no signal
        samp = length(x);
        y1(y1==0) = [];
        raw = y1;
        noiseinit=stdev;
        % store normalization constants
        M1(i)=mean(y1);
        S1(i)=std(y1);
        noiseinit=noiseinit/std(y1); 
        y1 = (y1 - mean(y1))/std(y1);
        raw = y1;
        %detrend data IMPORTANT PARAMETER - number (3rd input) controls how slow
        [m,par0] = detrenddataNEW(raw,x,DetrendParam);
        y1 = y1-m; %detrended y1
        % save new normalization constants
        M2(i)=mean(y1);
        S2(i)=std(y1);
        noiseinit=noiseinit/std(y1);
        y1 = (y1 - mean(y1))/std(y1); % mean normalised
        detrendData(i,:)= bring_to_size(y1',[1,numel(time)],NaN);
        % fit OU and OUoscillatory models
        [LLR, ll1,ll2,par1, par2,md] = fitDataLLR_2D(x,y1,noiseinit);
        par1TOT(i,:) = par1; % OU
        par2TOT(i,:) = par2; % OUosc
        LLRTOT(i,:) = LLR;
        LL1(i)=ll1;
        LL2(i)=ll2;
        dift(i).x=x;
        dfit(i).M=m*S1(i)+M1(i);
        dfit(i).m=md*S2(i)*S1(i)+M2(i);
        % in addition store these
        dfit(i).Raw=raw*S1(i)+M1(i);
        dfit(i).Detrended=y1*S1(i)*S2(i);
        modelFit(i,:)=bring_to_size(md',[1,numel(time)],NaN);
        showfigure(x,m,raw,y1,md,M1(i),S1(i),S2(i),LLR,par2,i);
        print(gcf(),[dirname,'/Cell',num2str(i)],'-dpng');
    end
    %% Step 4. Generate synthetic null scores
    repeats = round(200/length(par1TOT)); %to create synthetic timeseries - this can be increased or decreased
    [LLRM,LL1s,LL2s,par1Ms,par2Ms,synthOU]=getSynthLLRsimple(data,time,par1TOT,par0TOT,repeats);
    %% Step 5. Determine oscillators by FDR - first need to find pi_0
    % pi_0 compares shape of data with the synth OU to get prop of osc cells
    upper = max([LLRTOT;LLRM]);
    lower1 = min([LLRTOT;LLRM]);
    lower = upper - 0.90*(upper-lower1);
    range = linspace(lower,upper,20);
    for i= 1:length(range)
        cutoff = range(i);
        num = sum(LLRTOT<cutoff)/length(LLRTOT); % number of nonosc cells
        denom = sum(LLRM<cutoff)/length(LLRM); % number of nonosc cells in synthetic OU data
        piest(i) =  num/denom;
    end
    figure()
    subplot(1,2,1)
    plot(range,piest)
    ylim([0 1])
    hold off
    xlabel('\lambda')
    ylabel('\pi_0(\lambda)')
    xlim([lower1,upper])
    idx=find(piest==0);
    piest(idx)=[];
    range(idx)=[];
    %cubic spline regression
    xx = linspace(lower,upper,100);
    yy = spline(range,piest,xx);
    subplot(1,2,2),plot(xx,yy)
    % plot piGuess
    hold on
    plot(xx,yy,'color','r')
    ylim([0 1])
    hold off
    xlabel('\lambda')
    ylabel('\pi_0(\lambda)')
    xlim([lower1,upper])
    a = xlim();
    b = ylim();
    text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'a)'},...
        'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
    piGUESS1= min(1,min(yy));
    if piGUESS1<=0
        piGUESS1=1; % do not accept 0
    end
    %% Go through LLR scores calculating q values
    [LLRTOTs,I] = sort(LLRTOT);
    q1 = zeros(length(LLRTOT),1);
    for i = 1:length(LLRTOTs)
        Thresh = round(LLRTOTs(i));
        q1(i) = piGUESS1*(sum(LLRM>=Thresh)/length(LLRM))/(sum(LLRTOTs>=Thresh)/length(LLRTOTs));
    end
    
    cutoff = find(q1<q,1,'first')
    [w,l] = sort(I);
    Reorderedq = q1(l);
    Thresh=LLRTOTs(cutoff)
    try
        PassList = LLRTOT>Thresh;
    catch
        disp('No oscillators are found');
        PassList=zeros(size(LLRTOT));
    end
    %% plots LLR of data vs synthetic LLR
    subplot(2,1,1)
    m=max([LLRTOT(:);LLRM])+5;
    histogram(LLRTOT(:),15)
    t = title(fnames{k},'Interpreter','none');
    t.FontWeight = 'normal';
    xlim([0 m])
    subplot(2,1,2)
    histogram(LLRM,15)
    xlim([0 m])
    t = title('Synthetic bootstrap (non-osc)');
    t.FontWeight = 'normal';
    text(2,100,['Optimal LLR cutoff ', num2str(round(Thresh))]);
    print(gcf(),[dirname,'/LLR Summary'],'-dpng');
    %% Step 6. Summary data plots
    % analyse periods
    periods = 2*pi()./par2TOT(:,2);
    if ~isempty(Thresh)
        tmp=periods(PassList);
        figure,histogram(tmp,5)
        set(gca,'XLim',[0,10]);
        title(['Periods of passing cells ',num2str(round(mean(tmp)*100)/100),'h'])
        xlabel('Period (hours)')
        ylabel('Number of cells')
        print(gcf(),[dirname,'/PeriodPlot'],'-dpng');
        %% plots distribution of quality parameter of all cells
        Quality=par2TOT(:,2)./(2*pi*par2TOT(:,1));
        figure,histogram(Quality(PassList),5)
        title('Quality of oscillations')
        ylabel('Number of cells')
        print(gcf(),[dirname,'/QualityPlot'],'-dpng');
    end
    %% Step 7. Hilbert Analysis with Fold change detection
    hbt=HilbertFoldChange(dfit,dirname);
    %% extract selected stats from Hilbert analysis
    close all
    fcmax=[]; peak_count=[]; ptp_var=[];
    fc=[hbt(:).FoldChangeFit];
    for i=1:numel(hbt)
        vect=hbt(i).PairedIdx;
        peak_count(i)=size(vect,1);
        % peak to peak
        ppdist=(vect(2:end,1)-vect(1:end-1,1))*(time(2)-time(1));
        ptp_var(i)=var(ppdist);
        try
            fcmax(i)=max([hbt(i).FoldChangeFit]);
        catch
            fcmax(i)=NaN;
        end
    end
    figure,boxplot(fcmax),title('Maximum fold change')
    set(gca,'XTickLabel','');
    xlabel('All'), ylabel('Max Peak to Trough')
    print(gcf(),[dirname,'/Maximum fold change'],'-dpng');
    if ~isempty(Thresh)
        figure,boxplot(fcmax(PassList)),title('Maximum fold change in oscillators')
        set(gca,'XTickLabel','');
        xlabel('Oscillators'), ylabel('Max Peak to Trough');
        print(gcf(),[dirname,'/Maximum fold change osc'],'-dpng');
        list=1:numel(PassList);
        list(PassList)=[];
        figure,boxplot(fcmax(list)),title('Maximum fold change ')
        set(gca,'XTickLabel','');
        xlabel('Non-Oscillators'),ylabel('Max Peak to Trough');
        print(gcf(),[dirname,'/Maximum fold change nonosc'],'-dpng');
    end
    %% Step 8. Export data analysis and save simulation files
    s=size(data);
    % read in the cell information
    [~,str]=xlsread(fnames{k});
    str=str(2,coldata{k});
    str(rem)=[];
    % Column 1 track id
    savedata=zeros(s(2),10);
    header{1}='TrackId';
    % Column 2 processed cell id
    savedata(:,2)=1:s(2);
    header{2}='CellId';
    % Column 3 lengthscale
    savedata(:,3)=par2TOT(:,1);
    header{3}='Lengthscale(1/h)';
    % Column 4 period
    savedata(:,4)=2*pi./par2TOT(:,2);
    header{4}='Period(h)';
    % Column 5 LLR
    savedata(:,5)=LLRTOT;
    header{5}='LLR';
    % Column 6 quality
    savedata(:,6)=(par2TOT(:,2)./(2*pi*par2TOT(:,1)));
    header{6}='Quality';
    % Column 7 fod label 1=osc 0=nonosc
    savedata(:,7)=PassList;
    header{7}='FOD';
    % Column 8 max fold change
    savedata(:,8)=fcmax;
    header{8}='MaxFoldChange';
    % Column 9 number of peaks
    savedata(:,9)=peak_count;
    header{9}='NoPeaks';
    savedata(savedata>10^5)=Inf;
    % Column 10 peak to peak variance
    savedata(:,10)=ptp_var;
    header{10}='VarPeakToPeak';
    if ismac==0
        xlswrite([dirname,'/SummaryFile'],header,'Sheet1','A1:H1');
        xlswrite([dirname,'/SummaryFile'],savedata,'Sheet2','A2');
    else
        ExportWithHeader([dirname,'/SummaryFile_new.csv'],header,savedata(:,2:end),str);
    end
    %% save data
    save([dirname,'/',strtok(fnames{k},'.') ,'.mat'])
    %% rescale and export detrended and model fit
    detrendData=detrendData';
    detrendData=detrendData.*repmat(S1.*S2,size(detrendData,1),1);
    modelFit=modelFit';
    modelFit=modelFit.*repmat(S1.*S2,size(modelFit,1),1);
    if ismac==0
        xlswrite([dirname,'/ProcessedData.xlsx'],detrendData,'RawDetrendedRAW');
        xlswrite([dirname,'/ProcessedData.xlsx'],modelFit,'ModelFit');
    else
        csvwrite([dirname,'/ProcessedData_DetrendedData.csv'],detrendData);
        csvwrite([dirname,'/ProcessedData_ModelFit.csv'],modelFit);
    end
    %% export Hilbert stats over time
    FoldChangeMat=[];
    PeakMat=[];
    TroughMat=[];
    PhaseMap=[];
    for i=1:numel(hbt)
        vect=hbt(i).FoldChangeFit;
        FoldChangeMat=[FoldChangeMat,bring_to_size(vect',[100,1],NaN)];
        p=hbt(i).PairedIdx;
        PeakMat=[PeakMat, bring_to_size(p(:,1),[100,1],NaN)];
        TroughMat=[TroughMat,bring_to_size(p(:,2),[100,1],NaN)];
        ph=hbt(i).HilbertPhase;
        PhaseMap=[PhaseMap, bring_to_size(ph,[numel(time),1],NaN)];
    end
    if ismac==0
        xlswrite([dirname,'/Hilbert'],FoldChangeMat,'FoldChanges');
        xlswrite([dirname,'/Hilbert'],PeakMat,'HilbertPeaks');
        xlswrite([dirname,'/Hilbert'],TroughMat,'HilbertTroughs');
        xlswrite([dirname,'/Hilbert'],PhaseMap,'PhaseDiagram');
    else
        csvwrite([dirname,'/Hilbert_FoldChanges.csv'],FoldChangeMat);
        csvwrite([dirname,'/Hilbert_Peaks.csv'],PeakMat);
        csvwrite([dirname,'/Hilbert_Troughs.csv'],TroughMat);
        csvwrite([dirname,'/Hilbert_PhaseDiagram.csv'],PhaseMap);
    end
end
