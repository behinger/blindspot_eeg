function EEG = be_load_ICA(EEG,p,IDX)
%% Loads AMICA output specificed in p
% mv_load_ICA(EEG,p)
% Updates the Path
% Prints out all available ICA (if multiple)
% Loads the selected one, or if you provided an IDX loads the one on the
% IDX
% sample: mv_load_ICA(EEG,p,2) % adds the second ICA e.g. p.amica.path{2}
p = be_generate_paths(p);
if ~check_EEG(EEG.preprocess,'Ica')
    
    if length(p.amica.path) > 1
        for i = 1:length(p.amica.path)
            fprintf('%i : %s | %s \n',i,p.amica.date{i}, p.amica.path{i})
        end
        if nargin <3
            amica_index=input('Please choose an amica: ');
        else
            amica_index = IDX;
        end
    else
        amica_index = 1;
    end
    fprintf('loading AMICA: %s \n',p.amica.path{amica_index})
    addpath('/net/store/nbp/projects/EEG/blind_spot/amica')
    mod = loadmodout12(p.amica.path{amica_index});
    EEG.icaweights = mod.W;
    EEG.icasphere = mod.S;
    EEG.icawinv = [];EEG.icaact = [];EEG.icachansind = [];
    if isempty(findstr('Ica',EEG.setname))
        EEG.preprocess = [EEG.preprocess 'Ica'];
    end
    EEG.preprocessInfo.icaname = p.amica.path{amica_index};
    EEG.preprocessInfo.ICAloadDate = datestr(now);
    EEG.preprocessInfo.chosenICA = amica_index;
    EEG = eeg_checkset(EEG);
    %pop_saveset(EEG,'filename',EEG.setname,'filepath',p.path.set)
end
