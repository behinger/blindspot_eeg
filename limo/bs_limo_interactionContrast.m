function bs_limo_interactionContrast(LIMO)
%% runs contrast instead of anova
% flagIdx = 13;
% load(flags(flagIdx).group.LIMO)
% cd(flags(flagIdx).group.folder)
% figure,
% mat = bs_limo_plotDesignmatrix(LIMO);
cd([LIMO.dir '/groupLimo'])
load('LIMO.mat')

% pred1 = 5; %column
% pred2 = 7; %column
% intIdx = find(ismember(mat', (mat(:,pred1).*mat(:,pred2))','rows'));

tmp = load([LIMO.data.data_dir{1} '/LIMO.mat']);
% tmpIdx = find(strcmp(tmp.LIMO.design.XDesc,interactionLabel)); % find which predictor it is
% intIdx = sum(tmp.LIMO.design.nb_conditions)+tmpIdx-length(tmp.LIMO.design.nb_conditions); %column number in the designmatrix
endCondIdx = sum(tmp.LIMO.design.nb_conditions)+1;
if isempty(tmp.LIMO.design.nb_interactions) || (length(tmp.LIMO.design.nb_interactions) == 1 && tmp.LIMO.design.nb_interactions == 0)
    warning('no interaction found skipping')
else
    for k = 1:length(tmp.LIMO.design.nb_interactions)
        
        intIdx = endCondIdx+sum(tmp.LIMO.design.nb_interactions(1:(k-1))) : endCondIdx+sum(tmp.LIMO.design.nb_interactions(1:k)) - 1;
        if isempty(intIdx)
            error('interaction not found')
        end
        
        if length(intIdx) ~= 4
            continue
        end
        [~,pairs] = bs_limo_designMat_removeDouble(tmp.LIMO);
        dat1     = bs_limoGatherData('paths', LIMO.data.fullPath,'elec',[1:64],'predictor',intIdx);
        
        doublePreds = intIdx(ismember(intIdx,pairs));
        for w= 1:length(doublePreds)
            for x = 1:length(pairs)
                if any(doublePreds(w) == pairs(x,:))
                    newPred = pairs(x,mod(find(doublePreds(w) == pairs(x,:)),2)+1); % find and choose the other one
                    newRawDat = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',newPred);
                    newRawDat = permute(newRawDat,[1 2 4 3]); % just a trick to get the 3 dimension to be 1 for adding in the next line
                    dat1(:,:,doublePreds(w)==intIdx,:) = dat1(:,:,doublePreds(w)==intIdx,:)+newRawDat;
                end
            end
        end

        
        
        
        dat1 = permute(dat1,[1 2 4 3]);
        %     dat1 =
        dat2 = dat1(:,:,:,1) - dat1(:,:,:,2) - dat1(:,:,:,3) + dat1(:,:,:,4);
        LIMO.design.bootstrap = 1000;
        bs_limo_random_robust(1,dat2,100+k,LIMO.design.bootstrap,LIMO.design.tfce)
        
    end
end