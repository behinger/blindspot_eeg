% clear all
flags = be_check_folderstruct('bsLimo');
flagIdx = 19;

predList = 1
for subj = 1%1:10
    
    load(flags(flagIdx).subj(subj).LIMO)
    LIMO.data.fullPath = {flags(11).subj(subj).Betas};
    % LIMO
    LIMO.design.bootstrap = 600;
    
    bs_limo_display_results(1,['Condition_effect_' num2str(predList) '.mat'],flags(flagIdx).subj(subj).folder ,0.05,2,LIMO,...
        'title',sprintf('%s, Subj:%i',LIMO.design.XDesc{predList},subj),'addTopo','no','recalcMCC',0,'elecResort','no','caxis',[-5 5],'topoTime',210)
end

%% Interactions
predList = 1;
for subj = [4]
    
    load(flags(flagIdx).subj(subj).LIMO)
    LIMO.data.fullPath = {flags(flagIdx).subj(subj).Betas};
    % LIMO
    LIMO.design.bootstrap = 600;
    bs_limo_display_results(1,['Interaction_effect_' num2str(predList) '.mat'],flags(flagIdx).subj(subj).folder ,0.05,2,LIMO,...
        'title',sprintf('%s, Subj:%i',LIMO.design.XDesc{predList},subj),'addTopo','yes','recalcMCC',0,'elecResort','no','caxis',[-5 5])
end




%   bs_limo_display_results(1,['Rep_ANOVA_Factor_' num2str(predList) '.mat'],flags(flagIdx).group.folder ,0.05,2,LIMO,...
%         'title',sprintf('%s, Subj:%i','interact',subj),'addTopo','no','recalcMCC',0,'elecResort','no','caxis',[-150 150])
%% All in one plot
flagList =   {    19,       20,     };
predList    = {[1 2 3 4],   [1 2 3]};

for flagIdx = 1:length(flagList)
    for predictor = 1:length(predList{flagIdx})
%     predictor = 2;
%     flagIdx = 19;
%     figure,

    bs_limo_display_1stReduced(predList{flagIdx}(predictor),flags(flagList{flagIdx}))
    end
end