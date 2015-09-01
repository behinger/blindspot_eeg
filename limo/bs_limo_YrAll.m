for flagIdx = [55:57]
YrAll = [];
fprintf('\n\n *************\n loading flag:%d \n',flagIdx)
for k = 1:15
    fprintf('loading subject %d \n',k)
    load(flags(flagIdx).subj(k).Yr)
    load(flags(flagIdx).subj(k).LIMO)
    for l = 1:size(LIMO.design.X,2)
        if k == 1
            YrAll(:,:,l) = mean(Yr(:,:,LIMO.design.X(:,l)==1),3);
        else
        YrAll(:,:,l) = YrAll(:,:,l) + mean(Yr(:,:,LIMO.design.X(:,l)==1),3);
        end
    end
    
end
YrAll = YrAll/15;
save([flags(flagIdx).group.folder '/YrAll.mat'], 'YrAll')
end
