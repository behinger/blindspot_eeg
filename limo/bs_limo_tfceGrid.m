
for limoFlag = runForFlags(1:end)%[45 47 48]%39:43
    for pred = 1:5
%         continue
        cmd = sprintf('run(''~/blindspot.m'');eeglab;limo_eeg(-1);flags = be_check_folderstruct(''bsLimo'');bs_limo_display_resultsV2(%i,flags(%i).group.folder,''MCC'',5);',pred,limoFlag);
        nbp_grid_start_cmd('cmd',cmd,'name',sprintf('s_%i_f_%i',pred,limoFlag),'requ','exclusive=true,mem=1.5G','jobnum',randi(10000,1),'parallel',3,'queue','all')
        WaitSecs(15)
        while ~bs_qstat_check(20) %don't run more than 15 at the same time
            WaitSecs(30);
        end
    end
    for pred = 1:30

        cmd = sprintf('run(''~/blindspot.m'');eeglab;limo_eeg(-1);flags = be_check_folderstruct(''bsLimo'');bs_limo_display_resultsV2(%i,flags(%i).group.folder,''MCC'',5,''interaction'',''yes'');',pred+100,limoFlag);
        nbp_grid_start_cmd('cmd',cmd,'name',sprintf('s_%i_f_%i',pred+100,limoFlag),'requ','exclusive=true,mem=1.5G','jobnum',randi(100000,1),'parallel',3,'queue','all')
        WaitSecs(15)
        while ~bs_qstat_check(20) %don't run more than 15 at the same time
            WaitSecs(30);
        end
    end
end
