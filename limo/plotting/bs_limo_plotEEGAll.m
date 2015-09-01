%% Plot Main Effects

% flagList =   {    13,       14,        10,     15,       16,     8,   11 13};
% predList    = {[2 3 4],   [2 3 4],   [2 3],   [2 3],    [2 3], [2 3], [1,2], [1]};
% timeList    = {{[],[],[]},{[360],[],[]},{[],[]},{[],[]},{[],[]},{[],[]},{[],[210]},{[-20]}};
flags = be_check_folderstruct('bsLimo');
%%
% flagList =   {   22 23 24 25 26 27 29 30 31 32,33,34};
% flagList =   {   35 36 37 38};
flagList =   [42:43];
predList    = {[1 2 3],[1 2 3],[1 2 3],[1 2 3 4]};
elecList([39]).main = [22,7 ,32+5,32+18];
elecList([40]).main = [22,5+32 ,1,32+18];
elecList([41]).main = [22,7 ,1,32+18];
elecList([42]).main = [19+32,4,1,32+7];
elecList([43]).main = [22,7 ,1,32+7];


elecList(39).int = [30];
elecList(40).int = [7];
elecList(41).int = [30,7 ,1];
elecList(42).int = [5+32 32+7 ];
elecList(43).int = [5+32 32+7 ];
% timeList    = {{[],[]},{[],[],[],[]},{[],[],[]},{[],[],[]} {[],[],[],[]} {[],[],[],[]} {[],[],[],[]}};
plotFolder = '2014-04-22_pCorr'
plotFolder = '2014-04-22'
%%
%
alpha = 0.05;
alpha = 0.05/7;
badOnes = [];
g.print = 1;
for flagIdx = flagList
    for pred= 1:5
        try
%             bs_limo_display_resultsV2(pred,flags(flagIdx).group.folder,'MCC',5,'alpha',alpha)
             bs_limo_display_resultsV3(pred,flags(flagIdx).group.folder,elecList(flagIdx).main(pred),'MCC',5,'alpha',alpha,'style','')

            if g.print
                 be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/v3/%i/%i_%s',plotFolder,flagIdx,pred,flags(flagIdx).name),'prepare','plosBio','res','-r150','column',2)
%                  clf
                 
                 bs_limo_display_resultsV3(pred,flags(flagIdx).group.folder,elecList(flagIdx).main(pred),'MCC',5,'alpha',alpha,'style','topoPlotOnly')
                 be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/v3/%i/print/%i_%s_topo',plotFolder,flagIdx,pred,flags(flagIdx).name),'prepare','plosBio','res','-r300','column',2)
%                 bs_limo_display_resultsV3(pred,flags(flagIdx).group.folder,elecList(flagIdx).main(pred),'MCC',5,'alpha',alpha,'style','onlyCurves')
%                 be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/v3/%i/print/%i_%s_contour',plotFolder,flagIdx,pred,flags(flagIdx).name),'prepare','plosBio','res','-r150','column',2)
%                 
            end
            close all
        catch e
            try
                matlabpool('close')
            catch
            end
            badOnes{end+1} = {pred,flagIdx};
            warning(e.message)
%             e.stack(end)
        end
        
    end
end


%% Plot interactions
for flagIdx = flagList % for 1st Stim not defined
    for pred = 1:4
        try
            %             bs_limo_display_resultsV2(100+pred,flags(flagIdx).group.folder,'interaction','yes','MCC',5)
            
            try
                elecSelect = elecList(flagIdx).int(pred);
            catch
                elecSelect = 1;
            end
            bs_limo_display_resultsV3(100+pred,flags(flagIdx).group.folder,elecSelect,'MCC',5,'alpha',alpha,'interaction','yes')

            if g.print
                be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/v3/%i/%i_%s',plotFolder,flagIdx,100+pred,flags(flagIdx).name),'prepare','plosBio','res','-r150','column',2)
                
                bs_limo_display_resultsV3(100+pred,flags(flagIdx).group.folder,elecSelect,'MCC',5,'alpha',alpha,'interaction','yes','style','topoPlotOnly')
                 be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/v3/%i/print/%i_%s_topo',plotFolder,flagIdx,100+pred,flags(flagIdx).name),'prepare','plosBio','res','-r300','column',2)
            end
            close all
        catch e
            try
                matlabpool('close')
            catch
            end
            badOnes{end+1} = {pred,flagIdx};
            warning(e.message)
        end
        
        
        
    end
end
%     be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/v2/%i/%s_interaction',plotFolder,flagIdx,flags(flagIdx).name),'prepare','yes','res','-r100')
%     close(gcf)


%% Plot topoplots
for flagC = 1:length(flagList)
    flagIdx = flagList{flagC};
    
    for predC= 1:length(predList{flagC})
        pred = predList{flagC}(predC);
        bs_limo_topoplots(flagIdx,pred)
        be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/topoplots/%i/%i_%s_topoplot',plotFolder,flagIdx,pred,flags(flagIdx).name),'prepare','yes','res','-r100')
        
        close(gcf)
    end
    for k = 1:4
        try
            bs_limo_topoplots(flagIdx,100+k)
            be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/topoplots/%i/inter_%i_%s_topoplot',plotFolder,flagIdx,k,flags(flagIdx).name),'prepare','yes','res','-r100')
            
        catch
        end
    end
end
%% Plot Designmatrices

for flagIdx = flagList
%     flagIdx = flagList{flagC};
    bs_limo_plotDesignmatrix(flags(flagIdx).group.LIMO)
    
    be_print('file',sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/designMat/%i_designmat',plotFolder,flagIdx),'prepare','yes','res','-r100')
    close(gcf)
    %     end
end

%%

