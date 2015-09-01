if ~exist('flags','var')
    flags = be_check_folderstruct('bsLimo');
end
if ~exist('EEG','var')
    eeglab
    limo_eeg(0)
end

load(flags(37).group.LIMO)
% [~,L] = bs_limo_plotDesignmatrix(LIMO);

time = 300:400;
allDat = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',[5,6,23]);
allDat = squeeze(allDat(:,time,1,:) - allDat(:,time,2,:) + allDat(:,time,3,:));
allDat = mean(allDat,3);
%%

for k = 1:15
    fprintf('calc subj %i \n',k)
    subjDat = load(flags(37).subj(k).Yr);
    tmp = load(flags(37).subj(k).LIMO);
    subjDat.LIMO = tmp.LIMO;
    select = subjDat.LIMO.design.X(:,6)==1 & subjDat.LIMO.design.X(:,4)==1;
    tmpSubjDat = permute(subjDat.Yr(:,time,select),[3 1 2]);
   
%     rhoAll{k} = corr(allDat(:),tmpSubjDat(:,:)');
    rhoAll{k} = tmpSubjDat(:,:)*allDat(:);
end


%%
figure
jisubplot(3,5)

for k = 1:15
    nextplot
    tmp = rhoAll{k};
%     tmp = conv(tmp,ones(1,20)/20);
    plot(tmp/(64*100));
    axis([0 500 -1 1])
    hline(0)
    
end
% tmp = tmp/15;
% plot(tmp)
% myfun = @(x)std(x)/15
% shadedErrorBar(1:922,tmpAll,{@mean,myfun})