function [EEG reject] = be_ICA_mark(EEG,p,varargin)
%% be_ICA_mark(EEG,p,1,manualRejStruct)
% this function can be used to load the saved badComponents e.g.:
%  [~,reject] = be_ICA_mark(EEG,p,1)
% add new RejectStruct
%  be_ICA_mark(EEG,p,reject,1)
% manually clean
%  be_ICA_mark(EEG,p)
% add the rejections to EEG
%  EEG = be_ICA_mark(EEG,p,1)
if nargin > 2
    for l = 1:length(varargin)
        if isnumeric(varargin{l})
            silent = 1;
        elseif isstruct(varargin{l})
            newReject = varargin{l};
        else
            error('unkown varargin in be_ICA_mark')
        end
    end
else
    silent = 0;
end

if nargin <1
    error('not enough input arguments to be_ICA_mark')
end

EEG = eeg_checkset(EEG);

if exist(p.full.badComp,'file')==2
    tmpRej = load(p.full.badComp);
    fnRej = fieldnames(tmpRej);
    if isnumeric(tmpRej.(fnRej{1}))
        %we still have reject numbers
        reject.loaded = zeros(1,size(EEG.icawinv,2));
        reject.loaded(tmpRej.(fnRej{1})) = 1;
    else
        %we already have a reject struct
        reject =  tmpRej.(fnRej{1});
    end
end
if exist('newReject','var')
    fprintf('New RejectStruct loaded')
    reject = newReject;
    
end
if exist('reject','var') && isstruct(reject)
    allReject = any(cell2mat(struct2cell(reject)),1);

else
    fprintf('old file detected overwriting this session with nans \n')
    allReject = nan(1,size(EEG.icawinv,2));
end

tmp = dir(p.full.badComp);
fprintf('Reject Structure: \n',tmp.date)
fprintf('%i,',find(allReject==1)),fprintf('\n')

if ~silent
    askAppendOverwrite = input('Manualy clean? (y)/(n): ','s');
else
    askAppendOverwrite = 'silent';
end




switch askAppendOverwrite
    case 'y'
        EEG.reject.gcompreject = EEG.reject.gcompreject==1 | allReject == 1;
        EEG = pop_selectcomps_behinger(EEG,1:EEG.nbchan);
        uiwait;            fprintf('press any key to continue \n')            ,pause
        reject.manual = EEG.reject.gcompreject;
        resave = 1;
    case {'n','silent'}
        EEG.reject.gcompreject = allReject;
        resave = 0;
    otherwise
        error('User Canceled \n')
end
if exist('reject','var') && exist('newReject','var')
    fprintf('New Automatic Reject Struct found, overwriting \n')
    resave = 1;
end

if exist(p.full.badComp,'file')==2 && resave
    copyfile(p.full.badComp,[p.full.badComp '.bkp' datestr(now,'mm-dd-yyyy_HH-MM-SS')]);
    fprintf('Backup created \n')
end
if resave
    save(p.full.badComp,'reject');
    fprintf('Components Saved \n')
end
fn = fieldnames(reject);
for k = fn'
    EEG.reject.(k{:}) = reject.(k{:});
end
% EEG.reject.eTpruned = reject.eTpruned;
% EEG.reject.eTall = reject.eTall;
% EEG.reject.muscleRej = reject.muscleRej;
