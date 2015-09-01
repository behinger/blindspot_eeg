function [EEG] = be_epoch(EEG,p,varargin)
%% be_epoch(EEG,p,varargin)
% varargin can be:
% be_epoch(EEG,p,[triggerName],[epochTime],[bsl])
if nargin < 2
    error('wrong input number in be_epoch')
end
if nargin > 2 && ~isempty(varargin{1})
    trigger = varargin{1};
    if ~iscell(trigger)
        trigger ={trigger};
    end
else
   trigger = {'etFixOnset'};
end
if nargin > 2 && ~isempty(varargin{2})
    epochTime = varargin{2};
else
    epochTime = [-1 1.5];
end
if nargin > 3 && ~isempty(varargin{3})
    baseline = varargin{3};
else
    baseline = [-200 -50];
end
% TODO: Varargin --> RMbaseline, Epochtimes etc.
% if ~check_EEG(EEG.preprocess,'Epoch')
    EEG = pop_epoch( EEG, trigger, epochTime,  'epochinfo', 'yes','newname',EEG.preprocess);  
    EEG = pop_rmbase(EEG,baseline);
    EEG.preprocess = [EEG.preprocess 'Epoch'];
    EEG.preprocessInfo.epoch.trigger = trigger;
    EEG.preprocessInfo.epoch.epochTime = epochTime;
    EEG.preprocessInfo.epoch.bsl = baseline;
    EEG = eeg_checkset(EEG);
    %pop_saveset(EEG,'filename',EEG.preprocess,'filepath',p.path.set,'savemode','twofiles')
% end