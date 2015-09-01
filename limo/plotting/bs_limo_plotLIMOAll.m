%% Plot Main Effects

% flagList =   {    13,       14,        10,     15,       16,     8,   11 13};
% predList    = {[2 3 4],   [2 3 4],   [2 3],   [2 3],    [2 3], [2 3], [1,2], [1]};
% timeList    = {{[],[],[]},{[360],[],[]},{[],[]},{[],[]},{[],[]},{[],[]},{[],[210]},{[-20]}};
flagList =   {    19,       20,     21 };
predList    = {[1 2 3 4],   [1 2 3],   [1 2 3]};
timeList    = {{[],[],[],[]},{[],[],[]},{[],[],[]}};



for flagC = 1:length(flagList)
    flagIdx = flagList{flagC};

    for predC= 1:length(predList{flagC})
        pred = predList{flagC}(predC);
        load(flags(flagIdx).group.LIMO),
        bs_limo_plotTtest(LIMO,pred,'printPath',sprintf('/net/store/nbp/EEG/blind_spot/plots/27Jan14/2ndStim/%i',flagIdx),'time',timeList{flagC}{predC})
        close(gcf)
    end
end
%% Plot interactions ANOVA
flagList = {19};

for flagC = 1:length(flagList)
    flagIdx = flagList{flagC};
    
%     for pred= predList{flagC}
        load(flags(flagIdx).group.LIMO),
        
%         bs_limo_interactionANOVA(LIMO)
        bs_limo_interactionANOVAPlot(LIMO)
%         betterFilename = pred{1};
%         betterFilename(strfind(betterFilename,'/')) = '-';
        be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/27Jan14/2ndStim/interactions/%s_interaction',flags(flagIdx).name),'prepare','yes','res','-r100')
        close(gcf)
%     end
end

%% Plot interactions ONE SAMPLE TTEST
flagList = {13 15 16};
predList = {{'BS/Tempospatio Int','noTemp/noSpat Int', 'noTemp/Spat Int','Temp/noSpat Int','Temp/Spat Int'},...
            {'noTemp/noSpat Int', 'noTemp/Spat Int','Temp/noSpat Int','Temp/Spat Int'}, ...
            {'oBS/noTempospatio Int', 'oBS/Tempospatio Int','iBS/noTempospatio','iBS/Tempospatio'}, ...
            {'BS/Tempospatio Int'}};
for flagC = 1:length(flagList)
    flagIdx = flagList{flagC};
    if flagIdx < 14 % to plot only 
        continue
    end
    for pred= predList{flagC}
        load(flags(flagIdx).group.LIMO),
        bs_limo_interactionPlot(LIMO,pred{1})
        betterFilename = pred{1};
        betterFilename(strfind(betterFilename,'/')) = '-';
        be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/27Jan14/2ndStim/interactions/%s_interaction%s',flags(flagIdx).name,betterFilename),'prepare','yes','res','-r100')
        close(gcf)
    end
end

%% Plot Designmatrices

flagList =   [19 20 21]%[14 15 16];%[ 8 10 11 13 14 15 16]

for flagC = 1:length(flagList)
    flagIdx = flagList(flagC);
        bs_limo_plotDesignmatrix(flags(flagIdx).group.LIMO)
        
        be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/27Jan14/2ndStim/designMat/%i_designmat',flagIdx),'prepare','yes','res','-r100')
        close(gcf)
%     end
end


%%

