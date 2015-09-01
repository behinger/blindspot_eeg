function [] = bs_limo_ttestContrast(LIMO,k1,k2)
fullPath=LIMO.data.fullPath;
dat1 = bs_limoGatherData('paths', fullPath,'elec',[1:64],'predictor',k1);
dat2 = bs_limoGatherData('paths', fullPath,'elec',[1:64],'predictor',k2);


[~,LI] = bs_limo_plotDesignmatrix(LIMO);

[~,pairs] = bs_limo_designMat_removeDouble(LI);

dat1 = better_data(dat1,k1,pairs);
dat2 = better_data(dat2,k2,pairs);


bs_limo_random_robust(3,dat1,dat2,[k1 k2],LIMO.design.bootstrap,LIMO.design.tfce)

    function dat = better_data(dat,k,pairs)
        
        doublePreds = k(ismember(k,pairs));
        
        for w= 1:length(doublePreds)
            for x = 1:length(pairs)
                if any(doublePreds(w) == pairs(x,:))
                    newPred = pairs(x,mod(find(doublePreds(w) == pairs(x,:)),2)+1); % find and choose the other one
                    newRawDat = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',newPred);
                    dat = dat(:,:,:)+newRawDat;
                end
            end
        end
    end
end