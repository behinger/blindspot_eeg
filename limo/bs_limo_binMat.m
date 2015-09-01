
LIMO.dir = '/net/store/nbp/EEG/blind_spot/data/limo10-07_noInteraction';
[dirPath,dataPath,fullPath] = bs_limo_getFilepaths('filename','Condition_effect_3.mat','subjects',1:10,'path',LIMO.dir);


binMatRaw = zeros(64,400,10);
binMat = zeros(64,400,10);
for k = 1:length(fullPath)
    cPath = fullPath{k}; %current path
    load(cPath) 
    binMatTmp = zeros(64,400);
    binMatTmp(Condition_effect(:,:,2)<0.05) = 1;
%     binMatTmpRaw() = 1;
    binMat(:,:,k) = binMatTmp;
    binMatRaw(:,:,k) = Condition_effect(:,:,2);
end
%%
figure, imagesc(binMat(:,:,1))    
figure, imagesc(binMatRaw(:,:,1))    
binoinv([0.05 0.95],10,0.4013)
% bs_limo_display_results(1,[pairedFilename],pwd,0.05,2,LIMO,'interactive','yes','title','cond1Effect','addTopo','yes','caxis',[-20 20])
figure,imagesc(mean(binMat,3))
figure,imagesc(mean(binMat,3)>0.2)



figure, imagesc(mean(binMatRaw,3))
%%
for t = 1:400
    for e = 1:64
        histData = hist(squeeze(binMatRaw(e,t,:)),linspace(0,1,10));
        
        
        binUn(e,t) = 1-binocdf(sum(binMat(e,t,:)),10,0.05);
        chiUn(e,t) = 1-chi2cdf(sum((histData-1).^2/(1/10)),9);
        [~,ksUn(e,t)] = kstest2(histData,linspace(0,1,10));
    end
end
figure,imagesc(chiUn)
figure,imagesc(ksUn)
figure,imagesc(ksUn<0.01)
figure,imagesc((binUn))
figure,imagesc(log(binUn))
%%
binMatRaw(1,1,:)
figure, imagesc(binMatRaw(:,:,1))


%%
title('mean P over subjects')
title('thresholded mean P with >=70% of subjects individual p>=0.05')
title('mean over thresholded maps of subject individual p>=0.05')
title('ks test on uniformity of p-values')

title(['binomial probability to have x many significant subjects' 10 ' with 95% CI of binomial cumulative distribution with p = 0.05'])
title(['LOGARITHMIC binomial probability to have x many significant subjects' 10 ' with 95% CI of binomial cumulative distribution with p = 0.05'])