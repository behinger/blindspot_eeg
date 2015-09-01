function bs_limo_interactionANOVA(LIMO)
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
% tmpIdx = find(strcmp(tmp.LIMO.design.XDesc,interactionLabel)); % find which predictor it is
% intIdx = sum(tmp.LIMO.design.nb_conditions)+tmpIdx-length(tmp.LIMO.design.nb_conditions); %column number in the designmatrix
intIdx = sum(tmp.LIMO.design.nb_conditions)+1 : sum(tmp.LIMO.design.nb_conditions)+tmp.LIMO.design.nb_interactions;
if isempty(intIdx)
    error('interaction not found')
end
dat1 = bs_limoGatherData('paths', LIMO.data.fullPath,'elec',[1:64],'predictor',intIdx);

    dat1 = permute(dat1,[1 2 4 3]);
    LIMO.design.bootstrap = 1000;
   bs_limo_random_robust(6,dat1,ones(size(dat1,3),1),4,LIMO.design.bootstrap,LIMO.design.tfce)

