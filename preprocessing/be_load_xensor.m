function [EEG chanlocs] = be_load_xensor(EEG,p)
addpath('/net/store/projects/move/move_svn/code/ressources/')
[electrodes elec] = read_xensor(p.full.xensor,1,0);
chanlocs = [];

for ind = 6:size(elec.label,1) 
        chanlocs(end+1).labels = electrodes.labels{ind} ;
        chanlocs(end).X = elec.pnt(ind,1);
        chanlocs(end).Y = elec.pnt(ind,2);
        chanlocs(end).Z = elec.pnt(ind,3);
end
tmpFid = {'Nz','LPA','RPA'};
for k = 1:3
    chanlocs(end+1).labels =tmpFid{k};
    chanlocs(end).X = elec.pnt(k,1);
    chanlocs(end).Y = elec.pnt(k,2);
    chanlocs(end).Z = elec.pnt(k,3);
    chanlocs(end).type = 'FID';
end

EEG.chanlocs = chanlocs;
EEG=pop_chanedit(EEG, 'convert',{'cart2all'});
EEG = eeg_checkset(EEG);

EEG.preprocessInfo.xensor.arguments = {p.full.xensor,1,0};
EEG.preprocess = [EEG.preprocess 'Xensor'];
% chanlocs_out = EEG.chanlocs;