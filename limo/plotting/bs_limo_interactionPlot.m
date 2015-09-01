function bs_limo_interactionPlot(LIMO,interactionLabel)
% flagIdx = 13;
% load(flags(flagIdx).group.LIMO)
% cd(flags(flagIdx).group.folder)
% figure,
% mat = bs_limo_plotDesignmatrix(LIMO);
cd([LIMO.dir '/groupLimo'])
load('LIMO.mat')

% pred1 = 5; %column
% pred2 = 7; %column
% intIdx = find(ismember(mat', (mat(:,pred1).*mat(:,pred2))','rows'));
tmp = load([LIMO.data.data_dir{1} '/LIMO.mat']);
tmpIdx = find(strcmp(tmp.LIMO.design.XDesc,interactionLabel)); % find which predictor it is
intIdx = sum(tmp.LIMO.design.nb_conditions)+tmpIdx-length(tmp.LIMO.design.nb_conditions); %column number in the designmatrix

if isempty(intIdx) || intIdx == 0
    error('interaction not found')
end
dat1 = bs_limoGatherData('paths', LIMO.data.fullPath,'elec',[1:64],'predictor',intIdx);
if ~exist(['one_sample_ttest_parameter_' num2str(intIdx) '.mat'],'file')
limo_random_robust(1,dat1,[intIdx],LIMO.design.bootstrap,LIMO.design.tfce)
end



bs_limo_display_results(1,['one_sample_ttest_parameter_' num2str(intIdx) '.mat'],pwd,0.05,2,LIMO,...
        'title',sprintf('interaction %s',interactionLabel),'addTopo','no','recalcMCC',0,'elecResort','no','caxis',[-10 10])
    
    