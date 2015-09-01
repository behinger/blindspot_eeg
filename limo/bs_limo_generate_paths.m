function [paths] = bs_limo_generate_paths(folder)
if isstruct(folder) %in this case we want to update p!
    % the input file is in this case p!
    folder = folder.folder;
    fprintf('Updating the Folderstructure \n')
end
paths.folder = folder;


fileList = {'LIMO','Res','R2','Betas','Yhat','Yr'};
for k = fileList
    [tmpFileName tmpFolder] = bs_checkFileFolder(paths.folder,[k{:} '.mat']);
    if isempty(tmpFileName)
        paths.(k{:}) = [];
    else
        paths.(k{:}) = [tmpFileName];
    end
end
if ~isempty(tmpFolder)
    
    
    paths.folder = tmpFolder;
end


paths.project = 'bsLIMO';

function [file folder] = bs_checkFileFolder(folder,fileName)
file = [];
rawPathDir =dir(folder);
for k = 1:length(rawPathDir)
    if rawPathDir(k).isdir
        continue
    end
    tmpXf = rawPathDir(k).name;
    if any(strfind(tmpXf,fileName))
        file= [tmpXf];
    end
end
if isempty(file)
    file = [];
    folder = folder(1:end-1);%make A1 to A
    rawPathDir =dir(folder);
    for k = 1:length(rawPathDir)
        if rawPathDir(k).isdir
            continue
        end
        tmpXf = rawPathDir(k).name;
        if any(strfind(tmpXf,fileName))
            file= [tmpXf];
        end
    end
    
end
if isempty(file)
    folder = [];
else
    file  = [folder filesep file];
end