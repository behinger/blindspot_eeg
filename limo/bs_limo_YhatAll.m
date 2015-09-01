flags = be_check_folderstruct('bsLimo');
for flagIdx = [56:61]
YhatAll = [];
fprintf('\n\n *************\n loading flag:%d \n',flagIdx)
for k = 1:15
    fprintf('loading subject %d \n',k)
    load(flags(flagIdx).subj(k).Yhat)
    load(flags(flagIdx).subj(k).LIMO)
    for l = 1:size(LIMO.design.X,2)
        if k == 1
            YhatAll(:,:,l) = mean(Yhat(:,:,LIMO.design.X(:,l)==1),3);
        else
        YhatAll(:,:,l) = YhatAll(:,:,l) + mean(Yhat(:,:,LIMO.design.X(:,l)==1),3);
        end
    end
    
end
YhatAll = YhatAll/15;
save([flags(flagIdx).group.folder '/YhatAll.mat'], 'YhatAll')
end
