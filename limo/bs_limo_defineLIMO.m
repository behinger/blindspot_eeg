function LIMO = bs_limo_defineLIMO(varargin)

g = be_inputcheck(varargin, ...
    {'subject','string','','';
    'dataDir','string','',''});
if ischar(g) || isempty(g.dataDir)
    error(g)
end
LIMO.Level = 1;
LIMO.data.data_dir = g.dataDir;

LIMO.data.data = [g.subject 'merge.set'];

EEG = pop_loadset('filename', LIMO.data.data,'filepath',LIMO.data.data_dir,'loadmode','info');
LIMO.data.chanlocs = EEG.chanlocs;
LIMO.data.start = EEG.xmin;
LIMO.data.trim1 = EEG.xmin*1000;
LIMO.data.trim2 = EEG.xmax*1000;
LIMO.data.end = EEG.xmax;
LIMO.data.sampling_rate = EEG.srate;
% LIMO.data.Cat = 4;%Categorical variable(s)
% LIMO.data.Cat = [2 2 2 2];
LIMO.data.Cont = 0;
if exist([LIMO.data.data_dir '/../neighbour/neighbourMat_' g.subject '.mat'],'file') ~=0
    load([LIMO.data.data_dir '/../neighbour/neighbourMat_' g.subject]);
else
    EEG = pop_loadset('filename', LIMO.data.data,'filepath',LIMO.data.data_dir);
    [neighbours,neighbourMat] = limo_get_channeighbstructmat(EEG,55);
end
LIMO.data.neighbouring_matrix = neighbourMat;
neighbourDir = [LIMO.data.data_dir '/neighbour'];
save([LIMO.data.data_dir '/../neighbour/neighbourMat_' g.subject],'neighbourMat')

LIMO.design.zscore = 0;
LIMO.design.method = 'OLS';
LIMO.design.type_of_analysis = 'Mass-univariate';
LIMO.design.bootstrap = 0;
LIMO.design.tfce = 0;
LIMO.design.status = 'to do';


%%