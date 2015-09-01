function bs_limo_plotTtest(LIMO,pred,varargin)
g = be_inputcheck(varargin, ...
    {'time','integer',[],[];
    'printPath','string','','';
    'caxis','integer',[],[-20 20];
    });

statCorr = 2;


% flags = be_check_folderstruct('bsLimo');
% cd(flags(flagIdx).group.folder);
cd([LIMO.dir '/groupLimo'])
load('LIMO.mat')
% load(flags(flagIdx).group.LIMO);
addpath('/home/student/b/behinger/Documents/MATLAB/eeglab_dev/plugins/limo_v1.4/limo_cluster_functions')
predListAll = {'12','34','56','78',' 910','1112'};
plotTitle = {'Moving  Left vs. Right','Blindspot Off vs. On','temporal Contrast Off vs. On','spatial Contrast Off vs. On'};
% caxisList = {[-20 20] [-10 10] [-10 10] [-20 20]};

    predList =  predListAll{pred};
    pairedFilename = 'paired_samples_ttest_parameter_';
    designMatDesc = load([LIMO.data.data_dir{1} '/LIMO.mat']);
    
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
    [~,name] = fileparts(LIMO.dir);
    plotTitle = [name '  ' plotTitle];
    
    %     try
    bs_limo_display_results(1,[pairedFilename predList '.mat'],pwd,0.05,2,LIMO,...
        'interactive','yes','title',plotTitle,'addTopo','yes','topoTime',g.time,'caxis',g.caxis,'recalcMCC',0,'elecResort','no')
    if ~isempty(g.printPath)
%         set(gcf,'Position',[1 145 1920 950])
%         pause
        be_print('file',sprintf('%s/%s_mcc%i',g.printPath,plotTitle,statCorr),'prepare','yes','res','-r100')
    end
    
    
    
% end