function [LIMO pairs] = bs_limo_designMat_findDouble(LIMO)
%% Test for problems with the designmat:
if LIMO.Level == 2;
   warning('Level 2 detected, trying to load first LIMO dataset')
   load([LIMO.data.data_dir{1} '/LIMO.mat'])
end
d = LIMO.design;
LIMO.design.Xorg = LIMO.design.X;
X = (LIMO.design.X);
XDesc = bs_limo_designMat_betterDesc(LIMO.design.XDesc);

pairs = [];
for k = 1:size(X,2)
    
    for t = k+1:size(X,2)
       if all(X(:,k) == X(:,t))
           pairs = [pairs; [k t]];
%            fprintf('double row found: %d - and %d \n',k,t)
%            X(:,k) = [];
%            break
           
       end
    end
    
end
LIMO.design.X = X;
%% Which Interaction is based on which factors?
% 
% k = 21;
% X = unique(LIMO.design.X,'rows');
% XDesc = bs_limo_designMat_betterDesc(LIMO.design.XDesc);
% % ~,~,~,matchstring] = regexp(LIMO.design.XDesc{1},'([^,]*)')
% if k > sum(LIMO.design.nb_conditions) %% we have an interaction here
%     for mE1 = 1:sum(LIMO.design.nb_conditions) %mainEffect
%         % Two Way
%         for mE2 = mE1:sum(LIMO.design.nb_conditions) 
%             if all(X(:,k) == X(:,mE1).*X(:,mE2))
%                 fprintf('found for k = %d, mE1: %s, mE2: %s \n',k,XDesc{mE1},XDesc{mE2})
%             end
%         end
%     end
% end
%         