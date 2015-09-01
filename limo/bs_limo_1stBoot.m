subjList = {'A','B','C','D','E','G','H','J','K','L','M','N','T','U','V'}
flags= be_check_folderstruct('bsLimo')
%% Run the semi partial correlations on the grid
    flags =  be_check_folderstruct('bsLimo');
for subj =1:15
    for flagIdx =29 %22:25;
         while ~bs_qstat_check(15)
                WaitSecs(120);
         end
            fprintf('running %i in flag %i \n',subj,flagIdx)
            cmd = ['run ~/blindspot;eeglab;limo_eeg(0);cd ' flags(flagIdx).subj(subj).folder ';load(''LIMO.mat'');limo_semi_partial_coef(LIMO);'];
            %     bs_grid_ready
            nbp_grid_start_cmd('cmd',cmd,'requ','mem_free=4G,exclusive=true','name',sprintf('sPc_%i',subj),'jobnum',randi(100000,1),'out','/net/store/nbp/EEG/blind_spot/gridOutput')
            WaitSecs(5)
    end
end
%% Run the bootstraping of H0 on the grid
% pathToLIMO = ['/net/store/nbp/EEG/blind_spot/data/limo12-07_Interaction/' subjList{2} '/'];
%Done: Subj 1, all idx
% started 03 Feb s2:5 idx = 23
for subj =1:15
    
    flags =  be_check_folderstruct('bsLimo');
    for flagIdx =26 %22:25;
        cd(flags(flagIdx).subj(subj).folder)
        tic
        bs_limo_eeg_4_parallel([])
        toc
        addpath('/net/store/nbp/EEG/nbp_grid_script/')
        % pathToLIMO = ['/net/store/nbp/EEG/blind_spot/data/limoTest2/' subjList{2} '/'];
        stepsize = 1;
        for k = 1:stepsize:64
            while ~bs_qstat_check(25)
                WaitSecs(120);
            end
            WaitSecs(10)
            fprintf('running %i till %i',k,k+stepsize)
            cmd = ['run ~/blindspot;eeglab;cd ' flags(flagIdx).subj(subj).folder '; bs_limo_eeg_4_parallel([' num2str(k:k+stepsize-1) ']);'];
            %     bs_grid_ready
            nbp_grid_start_cmd('cmd',cmd,'requ','mem_free=4G,exclusive=true','name',sprintf('bs_limo_%i',k),'jobnum',randi(100000,1),'out','/net/store/nbp/EEG/blind_spot/gridOutput')
            
        end
        
    end
    
end
%% Rerun missingedit bs

flags =  be_check_folderstruct('bsLimo');
for subj =1:15
    flagIdx = 25;
    cd(flags(flagIdx).subj(subj).folder)
%     bs_limo_eeg_4_parallel([])
    addpath('/net/store/nbp/EEG/nbp_grid_script/')
    % pathToLIMO = ['/net/store/nbp/EEG/blind_spot/data/limoTest2/' subjList{2} '/'];
    for k = 1:stepsize:64
        if ~exist(sprintf('H0/modelTmp%i.mat',k),'file')
            fprintf('%i,%i',subj,k)
        while ~bs_qstat_check(100)
            WaitSecs(30);
        end
        fprintf('running %i till %i',k,k+stepsize)
        cmd = ['run ~/blindspot;eeglab;cd ' flags(flagIdx).subj(subj).folder '; bs_limo_eeg_4_parallel([' num2str(k:k+stepsize-1) ']);'];       
        nbp_grid_start_cmd('cmd',cmd,'requ','mem_free=4G,exclusive=true','name',sprintf('b_s%i_f%i_%i',subj,flagIdx,k),'jobnum',randi(100000,1),'out','/net/store/nbp/EEG/blind_spot/gridOutput')
        WaitSecs(90)
        end
    end
end
%% Delete all H0 Files
 
flags =  be_check_folderstruct('bsLimo');
for subj =1:15
    
    for flagIdx = 22:25;
        cd(flags(flagIdx).subj(subj).folder)
        allFold = dir('H0');
        for k = find(cellfun(@(x)any(strfind(x,'modelTmp')),{allFold.name}))
            delete(sprintf('H0/%s',allFold(k).name));
        end
    end
end

%% load the stuff and sort it
for subj =4:9
    flags =  be_check_folderstruct('bsLimo');
    flagIdx = 25;
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
    
    
    
end
%%

%%
for k = [1:15]
     while ~bs_qstat_check(2)
            WaitSecs(60);
        end
cmd = sprintf('run ~/blindspot;bs_limo_1stBootCombine(%i);',k);
nbp_grid_start_cmd('cmd',cmd,'requ','mem_free=12G,exclusive=true','name',sprintf('limo_combine_%i',k),'jobnum',k,'out','/net/store/nbp/EEG/blind_spot/gridOutput')
        
end