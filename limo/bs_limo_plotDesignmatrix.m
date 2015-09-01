function varargout = bs_limo_plotDesignmatrix(LIMO)
%% [[Limo.X],[limo]] = bs_limo_plotDesignmatrix(LIMO)
if ischar(LIMO)
    load(LIMO)
end
if LIMO.Level == 2;
   warning('Level 2 detected, trying to load first LIMO dataset')
   load([LIMO.data.data_dir{1} '/LIMO.mat'])
end
if nargout == 0
figure,
imagesc(LIMO.design.X)
colormap(gray)
idx = 0.5;
% set(gca,'XTick',(1.5:sum([length(LIMO.design.nb_conditions) length(LIMO.design.nb_interactions)]))*2)
% set(gca
XTicksLabels = [];%num2cell(get(gca,'XTicklabel'));
% XTicks = [];
for k = 1:length(LIMO.design.nb_conditions)
    idx = LIMO.design.nb_conditions(k) + idx;
    if k == 1
        XTicks = 1.5;
    else
        XTicks(end+1) = XTicks(end) + LIMO.design.nb_conditions(k);
    end
    try
    XTicksLabels{k} = LIMO.design.XDesc{k};
    catch
    end
    vline(idx,'r')
    
end
% length(LIMO.design.nb_conditions)
if ~isempty(LIMO.design.nb_interactions)&&sum(LIMO.design.nb_interactions)>0
   for k = 1:length(LIMO.design.nb_interactions)
       
        
      XTicks(end+1) = XTicks(end) + LIMO.design.nb_interactions(k);
      idx = LIMO.design.nb_interactions(k) + idx;
      vline(idx,'b')
        try
            XTicksLabels{end+1} = LIMO.design.XDesc{k+length(LIMO.design.nb_conditions)};
        catch
        end
   end
    
end

set(gca,'XTick',XTicks)
set(gca,'XTicklabel',XTicksLabels)
xticklabel_rotate([],45,[])
% endTick = sum([LIMO.design.nb_conditions(LIMO.design.nb_conditions==2) LIMO.design.nb_interactions]);
% xticklabel_rotate([1.5:2:endTick (1:sum(LIMO.design.nb_conditions==1))+endTick],45,XTicksLabels)
 

else

    varargout{1} = LIMO.design.X;
    varargout{2} = LIMO;
end
% mat = LIMO.design.X;
% LIMO.design
% set(findobj('Linewidth'),'Linewidth',2)