function [EEG] = be_import(p)
%% Imports NBP EEG data
% also updates the path before doing it.
% names the EEG.setname to the filename
% EEG = mv_import_cnt(p) % with p = mv_generate_paths(..)
if nargin < 1 || nargin >1
    help mv_import_cnt
    error('Not enough or too many input arguments %i, should be 1',nargin)
end
% p = be_generate_paths(p);

if strcmp(p.full.raw(end-2:end),'cnt')
    addpath('/net/store/projects/move/move_svn/code/ressources/anteepimport1.09')
    fprintf('ANT Detected, importing it. make sure to have the anteepimport1.09 path added to your library!(check the NBP wiki!) \n')
    EEG = pop_loadeep(p.full.raw,'triggerfile','on');
    EEG.preprocess = [p.filename(1:end-4)];
elseif strcmp(p.full.raw(end-2:end),'set') ||strcmp(p.full.raw(end-2:end),'fdt')
    fprintf('EEGlab Set detected \n')
    EEG=  pop_loadset(p.full.raw);
    EEG.preprocess = [p.filename(1:end-4)];
elseif strcmp(p.full.raw(end-3:end),'vhdr')
    fprintf('Brain Vision File detected \n')
    EEG = pop_loadbv(p.path.raw, p.filename);
    EEG.preprocess = [p.filename(1:end-5)];
end

EEG.preprocessInfo.import.Date = datestr(now);
tmp = dir(p.full.raw);
EEG.preprocessInfo.import.rawFileInfo = tmp;
EEG = eeg_checkset(EEG);