function [EEG] = be_clean_continuous(EEG,p,varargin)
%% Cleaning function
%Example Call
%  [EEG] = mv_clean_continuous(EEG,p)
% Checks whether file is already cleaned, asks user what to do if so
%
% [EEG] = mv_clean_continuous(EEG,p,1)
% Does not ask user and simply rejects the data
global rej
rej = [];
% If there is an additional Input, that is not a char, then go silently and
% simply use the available cleaning times. This is for grid use
if nargin == 3 && ~ischar(varargin{1})
    silent = varargin{1};
    cleaning_check = 0;
    
    %If an Char called 'Bene' is given, then we want to go into cleaning
    %check to register whether the data is cleaned good. "Qualitycheck"
elseif nargin == 3 && ischar(varargin{1}) && strcmpi(varargin{1},'bene')
    silent = 0;
    cleaning_check =1;
    fprintf('Welcome Bene to the cleaningCheck')
else
    %Else simply run the script
    silent = 0;
    cleaning_check = 0;
end

resave = 0;

% if ~check_EEG(EEG.preprocess,'Clean')
%IF not on grid, ask for cleanerName
if  silent == 0 && ~cleaning_check;
    cleanerName = input('Your name: ','s');
end
% Check if the file already exists
if exist(p.full.badCont,'file')==2
    %Load it
    tmpRej = load(p.full.badCont);
    fnRej = fieldnames(tmpRej);
    rej = tmpRej.(fnRej{1});
    if ~any(strcmp('sampRate',fnRej))%~exist('sampRate','var') 
        sampRate = p.defaultSampRate;
    else
        sampRate = tmpRej.sampRate;
    end
    fprintf('Rejections loaded \n')
    % Check whether an empty (temporary) file has been made
    if exist('tmp','var')
        warning('Loaded Rejections were Empty! MAYBE SOMEONE ELSE IS CURRENTLY CLEANING THE FILE')
    else %if not temporary, check Samplingrate
        
        if isempty(rej)
            warning('Rejection are empty please overwrite')
        else
            %                 if EEG.srate ~= sampRate
            %                     error('Cleaning Sampling Rate (%i) is different from EEG sampling Rate (%i)\n',sampRate,EEG.srate)
            %                 end
            if size(rej,2)-5 ~= EEG.nbchan
                rej = [rej(:,1:5) zeros(EEG.nbchan,size(rej,1))'];
                fprintf('Last time the Data had been cleaned with different number of channels, fixing it!\n ')
            end
        end
    end
    
    % If silent, the cleaningInput should be "use old"
    if silent
        cleaningInput='u';
        % if Cleaning Check, we want to see the cleaning Times and use
        % Append
    elseif cleaning_check
        cleaningInput='a';
    else
        %Else we ask what the user wants
        cleaningInput = input('Old cleaning times found: (o)verwrite, (u)se old cleaning, (a)ppend, (c)ancel y/n: ','s');
    end
    %If we overwrite, we simply ask the user to clean again
    if strcmp(cleaningInput,'o')
        eegplot_keyboard(EEG.data,'srate',EEG.srate,'winlength',8, ...
            'events',EEG.event,'wincolor',[1 0.5 0.5],'command','global rej,rej=TMPREJ',...
            'eloc_file',EEG.chanlocs);
        uiwait;
        resave = 1;
        %If we append, we use the old rejectiontimes for 'winrej'
    elseif strcmp(cleaningInput,'a')
        if isempty(rej)
            fprintf('The old Rejections were empty, (a)ppend does the same thing as (o)verwrite \n');
        else
            if size(rej,2)-5 ~= EEG.nbchan,
                error('This should not happen anymore, contact Bene rej: %i, cleaning: %i \n \n',size(rej,2)-5,EEG.nbchan);
                
            end
        end
        eegplot_keyboard(EEG.data,'srate',EEG.srate,'winlength',8, ...
            'events',EEG.event,'wincolor',[1 0.5 0.5],'command','global rej,rej=TMPREJ',...
            'eloc_file',EEG.chanlocs,'winrej',sort(rej)./(sampRate/EEG.srate));
        uiwait;
        
        if strcmp(input('Do you really want to append (y)es/(a)bort :','s'),'a')
            error('User Aborted')
        end
        resave = 1;
        
        % If we cancel, thwrow an error
    elseif strcmp(cleaningInput,'c') %
        error('User canceled the cleaning-action')
        % The only possible Input now can be 'u', for continuing the
        % cleaning
    elseif ~strcmp(cleaningInput,'u')
        error('User gave impossible input')
    end
    
    %If the file does not exist simply clean (same as 'o')
else
    tmp=[];
    cleaningInput = 'firstRun';
    if ~isdir(p.path.reject),          mkdir(p.path.reject);     end
    
    save(p.full.badCont,'tmp')
    eegplot_keyboard(EEG.data,'srate',EEG.srate,'winlength', 8,...
        'events',EEG.event,'wincolor',[1 0.5 0.5],'command','global rej,rej=TMPREJ'...
        ,'eloc_file',EEG.chanlocs);
    uiwait;
    resave = 1;
    % use this for temporary filter
    %eegplot(EEG_temp.data,'srate',EEG_temp.srate,'winlength',4,'events',EEG_temp.event,'wincolor',[1 0.5 0.5],'command','rej=TMPREJ');
end
if cleaning_check
    %deprevated
    %         %If the CleaningChecker said the cleaningtimes are OK, then save 1
    %         if exist('cleanerName','var')
    %             fprintf('Last Person Cleaning: %s \n',cleanerName)
    %         end
    %         cleaningTimesOK = strcmpi(input('Clean Enough(y/n)? ','s'),'y');
    %         if ~cleaningTimesOK
    %             cleaningTimesComment = input('Comment on why not OK: ','s');
    %         else
    %             cleaningTimesComment = [];
    %         end
    %
    %         load(p.full.badCont)
    %         fprintf('All Changes from this Session Removed \n')
    %         copyfile(p.full.badCont,[p.full.badCont '.bkp' datestr(now,'mm-dd-yyyy_HH-MM-SS')]);
    %         fprintf('Backup created \n')
    %         if ~exist('cleanerName','var')
    %             cleanerName = 'unknown';
    %         end
    %         save(p.full.badCont,'rej','sampRate','cleanerName','cleaningTimesOK','cleaningTimesComment');
    %         fprintf('Quit \n')
    %         return
else
    
    % save the times of the rejection
    if exist(p.full.badCont,'file')==2 && ~silent && ~strcmp(cleaningInput,'u') || resave % if we find the file, and we are not on the grid and we do not continue without cleaning anything, backup!
        copyfile(p.full.badCont,[p.full.badCont '.bkp' datestr(now,'mm-dd-yyyy_HH-MM-SS')]);
        fprintf('Backup created \n')
    end
    % We get the current samplingRate
    if exist('sampRate','var')
        sampRate_org = sampRate;
    else
        sampRate_org = EEG.srate;
    end
    sampRate = EEG.srate;
    %If not on grid and not just cleaning the Dataset, save it!
    if ~silent && ~strcmp(cleaningInput,'u') || resave
        if isempty(rej)
            error('rejection was empty! Did not save and aborted')
        end
        
        save(p.full.badCont,'rej','sampRate','cleanerName');
        fprintf('Rejections saved \n')
    end
    
    % Convert and reject the marked rejections
    tmprej = eegplot2event(rej, -1);
    tmprej(:,3) = tmprej(:,3) - 0.5*EEG.srate; %assuming cutoff -6dB is 0.5Hz!
    tmprej(:,4) = tmprej(:,4) + 0.5*EEG.srate; %assuming cutoff -6dB is 0.5Hz!
    tmprej(tmprej(:,3)<1,3) = 1;
    tmprej(tmprej(:,4)>EEG.pnts,4) = EEG.pnts;
    % Sorting the values
    [~,sIdx] = sort(tmprej(:,3));
    tmprej = tmprej(sIdx,:);
    
    throw_out = nan;
    while ~isempty(throw_out)
        throw_out = [];
        for i = 1:size(tmprej,1)-1
            if tmprej(i,4)>=tmprej(i+1,3)
                throw_out=[throw_out i+1];
                tmprej(i,3) = min(tmprej(i,3),tmprej(i+1,3));
                tmprej(i,4) = max(tmprej(i,4),tmprej(i+1,4));
            end
            
            
        end
        tmprej(throw_out,:) = [];
    end
    
    if ~strcmp(cleaningInput,'a')
        tmprej(:,3:4) = round(tmprej(:,3:4)/(sampRate_org/EEG.srate));
    end
    [EEG LASTCOM] = eeg_eegrej(EEG,tmprej(:,[3 4]));
    
    %Save it in the preprocess unit.
    EEG.preprocessInfo.tmprej = tmprej;
    EEG.preprocessInfo.cleanContDate = datestr(now);
    EEG.preprocess = [EEG.preprocess 'Clean'];
    %     end
end
