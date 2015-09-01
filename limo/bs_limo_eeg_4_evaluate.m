function bs_limo_eeg_4_evaluate(LIMO,modelTmp)
% update the files to be stored on the disk
nboot = 599;
H0_Betas = NaN(64, size(modelTmp{1}.Betas,1), size(LIMO.design.X,2), nboot);
H0_R2 = NaN(64, size(modelTmp{1}.Betas,1), 3, nboot); % stores R, F and p values for each boot

if LIMO.design.nb_conditions ~= 0
    tmp_H0_Conditions = NaN(64, size(modelTmp{1}.Betas,1), length(LIMO.design.nb_continuous), 2, nboot);
end

if LIMO.design.nb_interactions ~=0
    tmp_H0_Interaction_effect = NaN(64, size(modelTmp{1}.Betas,1),length(LIMO.design.nb_interactions), 2, nboot);
end

if LIMO.design.nb_continuous ~= 0
    tmp_H0_Covariates = NaN(64, size(modelTmp{1}.Betas,1), LIMO.design.nb_continuous, 2, nboot);
end

tic
for e = 1:64
    electrode = e;
    fprintf('working on: elec %i, \t time taken: %.2f \n',e,toc)
    model = modelTmp{e};
    H0_Betas(electrode,:,:,:) = model.Betas;
    for B = 1:599 % now loop because we use cells
        H0_R2(electrode,:,1,B) = model.R2{B};
        H0_R2(electrode,:,2,B) = model.F{B};
        H0_R2(electrode,:,3,B) = model.p{B};
        
        if prod(LIMO.design.nb_conditions) ~=0
            if length(LIMO.design.nb_conditions) == 1
                tmp_H0_Conditions(electrode,:,1,1,B) = model.conditions.F{B};
                tmp_H0_Conditions(electrode,:,1,2,B) = model.conditions.p{B};
            else
                for i=1:length(LIMO.design.nb_conditions)
                    tmp_H0_Conditions(electrode,:,i,1,B) = model.conditions.F{B}(i,:);
                    tmp_H0_Conditions(electrode,:,i,2,B) = model.conditions.p{B}(i,:);
                end
            end
        end
        
        if  ~isempty(LIMO.design.nb_interactions)%LIMO.design.fullfactorial == 1
            if length(LIMO.design.nb_interactions) == 1
                tmp_H0_Interaction_effect(electrode,:,1,1,B) = model.interactions.F{B};% XXX  changed the last 1,1,: to 1,1,B
                tmp_H0_Interaction_effect(electrode,:,1,2,B) = model.interactions.p{B};% XXX 
            else
                for i=1:length(LIMO.design.nb_interactions)
                    tmp_H0_Interaction_effect(electrode,:,i,1,B) = model.interactions.F{B}(i,:);% XXX 
                    tmp_H0_Interaction_effect(electrode,:,i,2,B) = model.interactions.p{B}(i,:); % XXX 
                end
            end
        end
        
        if LIMO.design.nb_continuous ~=0
            if LIMO.design.nb_continuous == 1
                tmp_H0_Covariates(electrode,:,1,1,B) = model.continuous.F{B};
                tmp_H0_Covariates(electrode,:,1,2,B) = model.continuous.p{B};
            else
                for i=1:LIMO.design.nb_continuous
                    tmp_H0_Covariates(electrode,:,i,1,B) = model.continuous.F{B}(i,:);
                    tmp_H0_Covariates(electrode,:,i,2,B) = model.continuous.p{B}(i,:);
                end
            end
        end
    end
    %                 if e == 1 % XXX BEHINGER
    %                     cd H0
    %                     save boot_table boot_table
    %                     save H0_Betas H0_Betas -v7.3
    %                     save H0_R2 H0_R2 -v7.3
    %                 end
end
%             close(h)
warning on;

% save data on the disk and clear out
cd H0
save H0_Betas H0_Betas -v7.3
save H0_R2 H0_R2 -v7.3


if prod(LIMO.design.nb_conditions) ~=0
    for i=1:length(LIMO.design.nb_conditions)
        name = sprintf('H0_Condition_effect_%g',i);
        H0_Condition_effect = squeeze(tmp_H0_Conditions(:,:,i,:,:));
        save(name,'H0_Condition_effect','-v7.3');
        clear H0_Condition_effect
    end
    clear tmp_H0_Conditions
end

if  ~isempty(LIMO.design.nb_interactions)%LIMO.design.fullfactorial == 1
    for i=1:length(LIMO.design.nb_interactions)
        name = sprintf('H0_Interaction_effect_%g',i);
        H0_Interaction_effect = squeeze(tmp_H0_Interaction_effect(:,:,i,:,:));
        save(name,'H0_Interaction_effect','-v7.3');
        clear H0_Interaction_effect
    end
    clear tmp_H0_Interaction_effect
end

if LIMO.design.nb_continuous ~=0
    for i=1:length(LIMO.design.nb_continuous)
        name = sprintf('H0_Covariate_effect_%g',i);
        H0_Covariate_effect = squeeze(tmp_H0_Covariates(:,:,i,:,:));
        save(name,'H0_Covariate_effect','-v7.3');
        clear H0_Covariate_effect
    end
    clear tmp_H0_Covariates
end

clear electrode model H0_R2; cd ..
disp(' ');