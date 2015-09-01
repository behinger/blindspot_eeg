function data = bs_limo_gatherGA(LIMO,varargin)
%%   bs_limoGatherData(varargin)
%{'paths','cell',[],[];
%    'elec','integer',[1:64],[1:64];
%    'predictor','integer',[1:9],1:9;
g = be_inputcheck(varargin,...
    {'elec','integer',[1:64],[1:64];
    'predictor','integer',[1:9],1:9;
    'whatToLoad','string',[],'';
    'select','integer',[],[];
    });
if ischar(g)
    error(g)
end

if ~iscell(LIMO.data.data_dir)
   LIMO.data.data_dir = { LIMO.dir};
end
for k = 1:length(LIMO.data.data_dir)
    g.paths{k} =  [LIMO.data.data_dir{k} filesep g.whatToLoad];
end
% g.paths{k} =  [g.paths{k} filesep g.whatToLoad '_' num2str(g.predictor)];
if isempty(g.select)
    GAname = sprintf('%s/groupLimo/%s_%s.mat', LIMO.dir,g.whatToLoad,num2str(g.predictor,'%1.i-'));
    if exist(GAname,'file')
        load(GAname);
        data = squeeze(data(g.elec,:,:,:));
        return
    end
end
%check whether it was already calculated


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
        if isempty(g.select)
            sel = any(X(:,g.predictor) == 1,2);
        else
            sel = all(X(:,g.select) == 1,2);
        end
        data(:,:,i) = mean(tmp(:,:,sel),3);
    else
        error('unrecognized')
    end
    
    
end
if isempty(g.select)
    save(GAname,'data')
end
data = squeeze(data(g.elec,:,:));

% if ~isempty(g.elec)
%     data = squeeze(data);
%         data(1,1:size(tmp_data,1),size(tmp_data,1:size(tmp_data,2)) = tmp_data;
% end
end
function X = get_designmat(fn)
p = fileparts(fn);
[X] = bs_limo_plotDesignmatrix([p filesep 'LIMO.mat']);


end