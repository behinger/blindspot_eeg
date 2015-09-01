function bs_limo_interactionANOVAPlot(LIMO)
cd([LIMO.dir '/groupLimo'])
load('LIMO.mat')


bs_limo_display_results(1,['Rep_ANOVA_Factor_1.mat'],pwd,0.05,2,LIMO,...
        'title',sprintf('interaction 1'),'addTopo','no','recalcMCC',0,'elecResort','no','caxis',[-30 30])
    
    