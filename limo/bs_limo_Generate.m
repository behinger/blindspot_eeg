function bs_limo_Generate(subjToRun,flagIdxToRun,runlocal)
global EEG
if isempty(EEG)
eeglab;
end
local_path = which('limo_eeg');
root = local_path(1:max(find(local_path == filesep))-1);
addpath([root filesep 'limo_cluster_functions'])
addpath([root filesep 'help']);

%%
if nargin < 3
    runlocal = 1;
end
flags = be_check_folderstruct('bsLimo');
%%
for flagIdx = flagIdxToRun%18%:21
    subjList = {'A','B','C','D','E','G','H','J','K','L','M','N','T','U','V'};
    date = '2014-04-01'
%     date = '19-02-2014'
    %     subjList = {'M','N'};
    for subj = subjToRun%[2:10]%1:length(subjList)
        if ~isempty(strfind(flags(flagIdx).name,'1stStim')) && ~isempty(strfind(flags(flagIdx).name,'long')) 
            data_merge = sprintf('/net/store/nbp/EEG/blind_spot/data/study/%s-1stStim_long',date)
            
        elseif ~isempty(strfind(flags(flagIdx).name,'1stStim')) &&  ~isempty(strfind(flags(flagIdx).name,'avg'))
            data_merge = sprintf('/net/store/nbp/EEG/blind_spot/data/study/%s-1stStim_avg',date)
            
        elseif ~isempty(strfind(flags(flagIdx).name,'1stStim'))
            data_merge = sprintf('/net/store/nbp/EEG/blind_spot/data/study/%s-1stStim',date)
                  
        elseif ~isempty(strfind(flags(flagIdx).name,'Saccade'))&&  ~isempty(strfind(flags(flagIdx).name,'avg'))
            date = '2014-06-16'
            data_merge = sprintf('/net/store/nbp/EEG/blind_spot/data/study/%s-saccade_avg',date)
            
        elseif ~isempty(strfind(flags(flagIdx).name,'Saccade'))
            date = '2014-06-16'
            data_merge = sprintf('/net/store/nbp/EEG/blind_spot/data/study/%s-saccade',date)
            
        elseif (~isempty(strfind(flags(flagIdx).name,'Go')) &&  ~isempty(strfind(flags(flagIdx).name,'avg')))
            date = '2014-05-28'
            data_merge = sprintf('/net/store/nbp/EEG/blind_spot/data/study/%s-go_avg',date)
            
        elseif ~isempty(strfind(flags(flagIdx).name,'long'))
            data_merge = sprintf('/net/store/nbp/EEG/blind_spot/data/study/%s-2ndStim_long',date)
            
        elseif ~isempty(strfind(flags(flagIdx).name,'avg'))
            data_merge = sprintf('/net/store/nbp/EEG/blind_spot/data/study/%s-2ndStim_avg',date)
            
        else
            data_merge = sprintf('/net/store/nbp/EEG/blind_spot/data/study/%s-2ndStim',date)
        end
        LIMO = bs_limo_defineLIMO('subject',subjList{subj},'dataDir',data_merge); % takes some time, the first time, because of channNeighbourstruct it is saved though;
        %     LIMO.dir = ['/net/store/nbp/EEG/blind_spot/data/limo10-07_noInteraction/' subjList{subj} '/'];
        LIMO.dir = flags(flagIdx).subj(subj).folder;
        %     LIMO.dir = ['/net/store/nbp/EEG/blind_spot/data/limo12-07_noInteraction/' subjList{subj} '/'];
        %     LIMO.dir = ['/net/store/nbp/EEG/blind_spot/data/limo01-07/' subjList{subj} '/'];
        if ~exist(LIMO.dir,'dir'),mkdir(LIMO.dir),end
        LIMO.design.fullfactorial = 1;
        if strfind(flags(flagIdx).name,'full')
            if any(strfind(flags(flagIdx).name,'1stStim')) || any(strfind(flags(flagIdx).name,'Saccade')) || any(strfind(flags(flagIdx).name,'Go'))
                opt.stim = 1;
            else
                opt.stim = 2;
            end
            if strfind(flags(flagIdx).name,'NoInt')
                LIMO.design.fullfactorial = 0;
                LIMO = bs_limo_designMat(LIMO,'full','yes','stim',opt.stim,'interaction','no');
                
            elseif strfind(flags(flagIdx).name,'allInt')
                LIMO.design.fullfactorial = 1;
                LIMO = bs_limo_designMat(LIMO,'full','yes','stim',opt.stim,'interaction','full');
            else
                LIMO.design.fullfactorial = 0;
                LIMO = bs_limo_designMat(LIMO,'full','yes','stim',opt.stim,'interaction','yes');
            end
            
        else
            if strfind(flags(flagIdx).name,'oBS')
                opt.bs = 'no';
            elseif strfind(flags(flagIdx).name,'wBS')
                opt.bs = 'yes';
            else        error('problem, neither oBS nor wBS');     end
            
            if strfind(flags(flagIdx).name,'oS')
                opt.ses = 'no';
            elseif strfind(flags(flagIdx).name,'wS')
                opt.ses = 'yes';
            else        error('problem, neither oS nor wS');     end
            
            if any(strfind(flags(flagIdx).name,'1stStim')) || any(strfind(flags(flagIdx).name,'Saccade'))
                opt.stim = 1;
            else
                opt.stim = 2;
            end
            
            if strfind(flags(flagIdx).name,'noInt')
                opt.interaction = 'no';%       LIMO.design.fullfactorial=0;
            elseif strfind(flags(flagIdx).name,'manInt')
                opt.interaction = 'man';
            elseif strfind(flags(flagIdx).name,'allInt')
                opt.interaction = 'full';
            end
            LIMO = bs_limo_designMat(LIMO,'bs',opt.bs,'session',opt.ses,'stim',opt.stim,'interaction',opt.interaction); %0 --> no int, 1 interactions
            
        end
        
        %interactions
         if strfind(flags(flagIdx).name,'IRLS')
             LIMO.design.method = 'IRLS';
         end
        cd(LIMO.dir)
        save('LIMO','LIMO')
        
        %do glm fit
        if runlocal
            limo_eeg(4)
        end
        load('LIMO')
        limo_semi_partial_coef(LIMO);
    end
end
% limo_random_select.m

%% Run the semi partial correlations on the grid


%%
% for flagIdx = 7:10;
if ~runlocal
    for flagIdx = flagIdxToRun
        addpath('/net/store/nbp/EEG/nbp_grid_script/')
        
        for k =subjToRun
            while ~bs_qstat_check(10)
                WaitSecs(60);
            end
            cmd = ['run ~/blindspot;eeglab;' 10 'flags = be_check_folderstruct(''bsLimo'');'];
            cmd = [cmd '' sprintf('flagIdx=%i;load(flags(flagIdx).subj(%i).LIMO);cd(LIMO.dir);limo_eeg(4);',flagIdx,k) ''];
            
            nbp_grid_start_cmd('cmd',cmd,'name',sprintf('glmFit%i',k),'jobnum',k+10*flagIdx,'requ','exclusive=true','out','/net/store/nbp/EEG/blind_spot/gridOutput')
        end
    end
    % end
end