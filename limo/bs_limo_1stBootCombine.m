function bs_limo_1stBootCombine(subj)
eeglab
local_path = which('limo_eeg');
root = local_path(1:max(find(local_path == filesep))-1);
addpath([root filesep 'limo_cluster_functions'])
addpath([root filesep 'help'])


    flags =  be_check_folderstruct('bsLimo');
    flagIdx = 20;
    cd(flags(flagIdx).subj(subj).folder)
    
    allModel = {};
    allFold = dir('H0');
    for k = find(cellfun(@(x)any(strfind(x,'modelTmp')),{allFold.name}))
        
        tmpModel = load(sprintf('H0/%s',allFold(k).name));
        if ~isnan(str2double(allFold(k).name(end-5)))
            l = str2double(allFold(k).name(end-5:end-4));
        else
            l = str2double(allFold(k).name(end-4));
        end
        fprintf('loading %i-%i \n',l,length(tmpModel.modelTmp))
        
        allModel(l:length(tmpModel.modelTmp)) = tmpModel.modelTmp(l:end);
        %         delete(sprintf('H0/%s',allFold(k).name));
    end
    
    load(flags(flagIdx).subj(subj).LIMO)
    bs_limo_eeg_4_evaluate(LIMO,allModel)%
