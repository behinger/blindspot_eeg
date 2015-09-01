function [outE1 outE2] = bs_limo_designMat_getDesc(LIMO,k)
if LIMO.Level == 2;
    warning('Level 2 detected, trying to load first LIMO dataset')
    load([LIMO.data.data_dir{1} '/LIMO.mat'])
end
% if k > 100
%     %Interaction
%             k = sum(LIMO.design.nb_conditions)+(k-100)*4;
%
% end

X = unique(LIMO.design.X,'rows');
XDesc = bs_limo_designMat_betterDesc(LIMO.design.XDesc);
% ~,~,~,matchstring] = regexp(LIMO.design.XDesc{1},'([^,]*)')
outK1 = [];
intCumsum = [0 cumsum(LIMO.design.nb_interactions)];
kMain = sum(LIMO.design.nb_conditions);
idx = find((k - cumsum([LIMO.design.nb_conditions LIMO.design.nb_interactions]))<0,1,'first')-  length(LIMO.design.nb_conditions);
if k > sum(LIMO.design.nb_conditions) %% we have an interaction here
    for kI = intCumsum(idx)+1:intCumsum(idx+1)
        outE1 = [];
        outE2 = [];
        
        for mE1 = 1:sum(LIMO.design.nb_conditions) %mainEffect
            % Two Way
            for mE2 = mE1+1:sum(LIMO.design.nb_conditions)
                if LIMO.design.nb_interactions(idx) == 4
                    if all(X(:,kI+kMain) == X(:,mE1).*X(:,mE2))
                        fprintf('found for k = %d, mE1: %s, mE2: %s \n',k,XDesc{mE1},XDesc{mE2})
                        outE1 = [outE1 {XDesc{mE1} ',' XDesc{mE2}}];
                        %                 outE2 = [outE2 XDesc(mE2)];
                    end
                elseif LIMO.design.nb_interactions(idx) == 8
                    for mE3 = mE2+1:sum(LIMO.design.nb_conditions)
                        if all(X(:,kI+kMain) == X(:,mE1).*X(:,mE2).*X(:,mE3))
                            fprintf('found for k = %d, mE1: %s, mE2: %s \n',k,XDesc{mE1},XDesc{mE2})
                            outE1 = [outE1 {XDesc{mE1} ',' XDesc{mE2} ',' XDesc{mE3}}];
                            %                 outE2 = [outE2 XDesc(mE2)];
                        end
                    end
                end
            end
            
        end
        outK1{end+1} = ['||' unique(outE1)];
        
    end
end
% outE1 = unique([outK1{:}]);
outE1 = [outK1{:}];
outE2 = [];