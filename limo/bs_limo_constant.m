load([flags(55).group.LIMO])
[~,LI] = bs_limo_plotDesignmatrix(LIMO);

pred = size(LI.design.X,2);

constant = bs_limo_gatherGA(LIMO,'predictor',pred,'whatToLoad','Yr.mat');

%%

topoDiff = mean(constant,3);

timevec = round((LIMO.data.start:1/LIMO.data.sampling_rate:LIMO.data.end)*1000);
% topoDiff = rawData(:,:,1) - rawData(:,:,2);
figure
scale = max([abs(min(topoDiff(:))) abs(max(topoDiff(:)))]);
topoTimes = [-300:25:400];

jisubplot(5,8,0,'landscape',[0.1 0.2])
for k = 1:length(topoTimes)-1
    nextplot('byrow')
    fromT = find(timevec<=topoTimes(k),1,'last');
    toT = find(timevec<topoTimes(k+1),1,'last');
    tmpDat= mean(topoDiff(:,fromT:toT),2);
    topoplot(tmpDat,LIMO.data.chanlocs,'maplimits',[-scale scale],'nosedir','+Y');
    title(sprintf('%.0f <= %.0f',timevec(fromT),timevec(toT)))
    
    
end

cbar
figure
plot(timevec,topoDiff','k')
%%
figure,jisubplot(5,3)
for k = 1:15
    nextplot
    plot(timevec,constant(:,:,k)','k')
    axis([-300 500 -15 15])
end
%%
[m,i] = max(max(topoDiff,[],1) - min(topoDiff,[],1));

fprintf('Maximal Amplitude %.2f at  %i ms\n',m,timevec(i))