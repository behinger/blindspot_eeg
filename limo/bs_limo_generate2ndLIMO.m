function LIMO = bs_limo_generate2ndLIMO(destPath)


chanLoad = load('/net/store/nbp/EEG/blind_spot/data/channel_loc_default.mat');
% chanLoad = load(['/net/store/nbp/EEG/blind_spot/data/study/neighbour/neighbourMat_A']);
% exp_chanloc = load(['/net/store/nbp/EEG/blind_spot/data/expected_chanlocs']);

LIMO = [];
% LIMO.dir = '/net/store/nbp/EEG/blind_spot/data/limo12-07_Interaction';
LIMO.dir =   fileparts(destPath); %only the main path

[LIMO.data.data_dir,LIMO.data.data,LIMO.data.fullPath] = bs_limo_getFilepaths('filename','Betas.mat','subjects',1:15,'path',LIMO.dir);
for k = 1:5;fprintf('***************\n');end
fprintf('15 Subjects selected \n')
for k = 1:5;fprintf('***************\n');end
tmp = load([LIMO.data.data_dir{1} '/LIMO.mat']);
LIMO.Level = 2;
LIMO.data.chanlocs = fieldtripchan2eeglab(chanLoad.elec);
LIMO.data.neighbouring_matrix = chanLoad.elec.channeighbstructmat;
LIMO.data.sampling_rate = tmp.LIMO.data.sampling_rate;
LIMO.data.trim1= tmp.LIMO.data.trim1;
LIMO.data.start= tmp.LIMO.data.start;
LIMO.data.trim2= tmp.LIMO.data.trim2;
LIMO.data.end=tmp.LIMO.data.end;
LIMO.design.bootstrap = 1000;
LIMO.design.tfce = 0;
LIMO.design.analyseType = 'Full scalp analysis';
LIMO.design.electrode = [];
LIMO.design.x = [];
LIMO.design.name = 'paired t-test all elecs';
if ~exist(destPath,'dir'),mkdir(destPath);end
cd(destPath)
save(['LIMO.mat'],'LIMO')
fullPath = LIMO.data.fullPath;
save('fullpath.mat','fullPath')