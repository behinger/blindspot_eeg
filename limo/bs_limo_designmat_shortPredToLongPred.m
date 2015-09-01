function [predName] = bs_limo_designmat_shortPredToLongPred(LI,predOrg)
if predOrg>100
        
        tmpPred = sum(LI.design.nb_conditions)+sum(LI.design.nb_interactions(1:predOrg-100))-1;
        [E1,E2] = bs_limo_designMat_getDesc(LI,tmpPred);
        if isempty(E1)
            predName = 'unknown';
            return
        end
        if isempty(E2)
            predName = [E1{:}];
        else
        predName = [E1{:} E2{:}];
        end
    else
        tmpPred = predOrg*2;
        predName = [LI.design.XDesc{predOrg}];
 end
 
