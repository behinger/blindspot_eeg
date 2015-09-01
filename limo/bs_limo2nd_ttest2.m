%% Second Level
% [Names,Paths,Files] = limo_get_files
flagIdx = 18;
LIMO = bs_limo_generate2ndLIMO(flags(flagIdx).group.folder);

%%

for k = [1 3 5 7 9 11;2 4 6 8 10 12]
    dat1 = bs_limoGatherData('paths', LIMO.data.fullPath,'elec',[1:64],'predictor',[k(1)]);
    dat2 = bs_limoGatherData('paths', LIMO.data.fullPath,'elec',[1:64],'predictor',[k(2)]);
    tic
    bs_limo_random_robust(7,dat1,dat2,[k(1) k(2)],LIMO.design.bootstrap,LIMO.design.tfce)
    toc
%     limo_random_robuest(1,dat1,[18],LIMO.design.bootstrap,LIMO.design.tfce)

end
%% Grid version T-Tests
addpath('/net/store/nbp/EEG/nbp_grid_script/')
for flagIdx =18%19:21;
    
    for k = [1 3 5 7; 2 4 6 8]
        while ~bs_qstat_check(16)
            WaitSecs(120);
        end

        %     if k(1) == 1 &&load(flags(13).group.LIMO) flagIdx == 3;
        %         break
        %     end
        cmd = sprintf(['flags = be_check_folderstruct(''bsLimo'');fIdx = %i;cd(flags(fIdx).group.folder);load(flags(fIdx).group.LIMO);k1=%i;k2=%i;' 10 ...
            'fullPath=LIMO.data.fullPath; dat1 = bs_limoGatherData(''paths'', fullPath,''elec'',[1:64],''predictor'',k1);' 10 ...
            'dat2 = bs_limoGatherData(''paths'', fullPath,''elec'',[1:64],''predictor'',k2);' 10 ...
            'bs_limo_random_robust(7,dat1,dat2,[k1 k2],LIMO.design.bootstrap,LIMO.design.tfce)'],flagIdx,k(1),k(2));
        
        nbp_grid_start_cmd('cmd',['run ~/blindspot;eeglab;' cmd],'name','bs_l_ttest','requ','mem_free=3G,exclusive=true','jobnum',k(1)+flagIdx*8,'out','/net/store/nbp/EEG/blind_spot/gridOutput')
        %'mem_free=14G,exclusive=true'
    end
    
end
%% Grid Version Interaction ANOVAs
for flagIdx =19%:21;
cmd = sprintf(['flags = be_check_folderstruct(''bsLimo'');fIdx = %i;cd(flags(fIdx).group.folder);load(flags(fIdx).group.LIMO);bs_limo_interactionANOVA(LIMO)'],flagIdx);
        
nbp_grid_start_cmd('cmd',['run ~/blindspot;eeglab;' cmd],'name','bs_anova_interaction','requ','mem_free=3G,exclusive=true','jobnum',20+flagIdx*18,'out','/net/store/nbp/EEG/blind_spot/gridOutput')
% bs_limo_interactionANOVA(LIMO)


end
% bs_limo_display_results(1,'one_sample_ttest_parameter_18.mat','/net/store/nbp/EEG/blind_spot/data/limo/limo19_09wBSoS/groupLimo/',0.15,2,LIMO,'recalcMCC',0)
%%

flagIdx = 18;
int = 'no'; 
int = 'yes';
statCorr = 2;
print = 0;
topoTimeList = [-80 366 300 240];
% topoTimeList = [-80 -18 160 240];
% 365


% flags = be_check_folderstruct('bsLimo');
cd(flags(flagIdx).group.folder);
load(flags(flagIdx).group.LIMO);
addpath('/home/student/b/behinger/Documents/MATLAB/eeglab_dev/plugins/limo_v1.4/limo_cluster_functions')
predListAll = {'12','34','56','78',' 910','1112'};
plotTitle = {'Moving  Left vs. Right','Blindspot Off vs. On','temporal Contrast Off vs. On','spatial Contrast Off vs. On'}
caxisList = {[-20 20] [-10 10] [-10 10] [-20 20]};
for k = 3
    predList =  predListAll{k};
    pairedFilename = 'paired_samples_ttest_parameter_';
    designMatDesc = load(flags(flagIdx).subj(1).LIMO);
    if isfield(designMatDesc.LIMO.design,'XDesc')
        plotTitle = designMatDesc.LIMO.design.XDesc;
    end
    switch predList
        case '12'
            plotTitle = plotTitle{1};
        case '34'
            plotTitle = plotTitle{2};
        case '56'
            plotTitle = plotTitle{3};
        case '78'
            plotTitle = plotTitle{4};
    end
    plotTitle = [flags(flagIdx).name '  ' plotTitle];
    if print == 1;int = 'no'; end
    %     try
    bs_limo_display_results(1,[pairedFilename predList '.mat'],pwd,0.05,statCorr,LIMO,...
        'interactive',int,'title',plotTitle,'addTopo','yes','topoTime',topoTimeList(k),'caxis',caxisList{k},'recalcMCC',0,'elecResort','no')
    if print
        set(gcf,'Position',[1 145 1920 950])
        pause
        be_print('file',sprintf('/work/behinger/Dropbox/Masterarbeit/Plots/eeg/2ndLevel/%s_time%i_mcc%i',plotTitle,topoTimeList(k),statCorr),'prepare','yes','res','-r100')
    end
    
    
    
end