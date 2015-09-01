function data = bs_limoGatherData(varargin)
%%   bs_limoGatherData(varargin)
%{'paths','cell',[],[];
%    'elec','integer',[1:64],[1:64];
%    'predictor','integer',[1:9],1:9;
g = be_inputcheck(varargin,...
    {'paths','cell',[],[];
    'elec','integer',[1:64],[1:64];
    'predictor','integer',[1:9],1:9;
    'whatToLoad','string',[],'';
    'overwriteFilename','string',[],'';
    });
if ischar(g)
    error(g)
end
if ~isempty(g.whatToLoad) || ~isempty(g.overwriteFilename)
    if length(g.predictor) > 1
       error('having multiple predictors and custom data is not supported') 
    end
    for k = 1:length(g.paths)
        if ~isempty(g.overwriteFilename)
          g.paths{k} =  [g.paths{k} filesep g.overwriteFilename];
        else
          g.paths{k} =  [g.paths{k} filesep g.whatToLoad '_' num2str(g.predictor)];
        end
    end
end
% [first_frame,last_frame,subj_chanlocs,channeighbstructmat] = match_frames(LIMO2lvl.data.data_dir);


for i=1:size(g.paths,2) % for each subject
    tmp = load([g.paths{i}]);
    fnTmp = fieldnames(tmp);
    tmp = tmp.(fnTmp{1});
    
    
    %         data(:,:,g.predictor,i) = tmp(g.predictor;
    if findstr(g.paths{i},'Betas')
        data(:,:,1:length(g.predictor),i) = tmp(g.elec,:,g.predictor);
    elseif findstr(g.paths{i},'R2')
        data(:,:,i) = tmp(g.elec,:,1);
    elseif  findstr(g.paths{i},'semi')
         data(:,:,i) = tmp(g.elec,:,1);
     elseif  any(findstr(g.paths{i},'Yr')) || any(findstr(g.paths{i},'Yhat'))
         X = get_designmat(g.paths{i});
         sel = any(X(:,g.predictor) == 1,2);
         data(:,:,i) = mean(tmp(g.elec,:,sel),3);
    else
        error('unrecognized')
    end
    
    
    clear tmp
end

if ~isempty(g.elec)
    data = squeeze(data);
    %         data(1,1:size(tmp_data,1),size(tmp_data,1:size(tmp_data,2)) = tmp_data;
end
end
function X = get_designmat(fn)
p = fileparts(fn);
[X] = bs_limo_plotDesignmatrix([p filesep 'LIMO.mat']);


end