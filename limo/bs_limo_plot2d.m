
flags = be_check_folderstruct('bsLimo');
subj = 4
%%
for flagIdx = 3:6
    cd(flags(flagIdx).subj(subj).folder);
    load(flags(flagIdx).subj(subj).LIMO);
    fileName = 'R2';
    % figure
    bs_limo_display_results(1,[fileName '.mat'],pwd,0.05,1,LIMO,'interactive','no','title',[flags(flagIdx).name ' ' fileName],'addTopo','no','elecResort','yes','caxis',[-20 20])
end