%% Plot Main Effects

% flagList =   {    13,       14,        10,     15,       16,     8,   11 13};
% predList    = {[2 3 4],   [2 3 4],   [2 3],   [2 3],    [2 3], [2 3], [1,2], [1]};
% timeList    = {{[],[],[]},{[360],[],[]},{[],[]},{[],[]},{[],[]},{[],[]},{[],[210]},{[-20]}};
flags = be_check_folderstruct('bsLimo');
%%
% flagList =   {   22 23 24 25 26 27 29 30 31 32,33,34};
% flagList =   {   35 36 37 38};
% predList    = {[1 2 3],[1 2 3],[1 2 3],[1 2 3 4]};
% elecList([39]).main = [22,7 ,32+5,32+18];
% elecList([40]).main = [22,5+32 ,1,32+18];
% elecList([41]).main = [22,7 ,1,32+18];
% elecList([42]).main = [19+32,4,1,32+7];
% elecList([43]).main = [22,7 ,1,32+7];
%
%
% elecList(39).int = [30];
% elecList(40).int = [7];
% elecList(41).int = [30,7 ,1];
% elecList(42).int = [5+32 32+7 ];
% elecList(43).int = [5+32 32+7 ];


%%
flagList =   [44];

g.print = 0;
plotFolder = '2014-04-29';
for k = 2
    
    if k == 1
        alpha = 0.05;
    else
        alpha = 0.05/11;
    end
    
    filePath = sprintf('/net/store/nbp/EEG/blind_spot/plots/%s/alpha_%.4f/%%s/%%s_%%i',plotFolder,alpha);
    %
    % alpha = 0.05;
    
    badOnes = [];
    
    for flagIdx = flagList
        continue
        for pred= 1:5
            try
                %             bs_limo_display_resultsV2(pred,flags(flagIdx).group.folder,'MCC',5,'alpha',alpha)
                bs_limo_display_resultsV4(pred,flags(flagIdx).group.folder,[],'MCC',5,'alpha',alpha,'style','')
                load(flags(flagIdx).subj(1).LIMO);
                thisFileName = sprintf(filePath ,LIMO.design.XDesc{pred},flags(flagIdx).name,pred);
                if g.print
                    be_print('file',thisFileName,'prepare','plosBio','res','-r150','column',2)
                    
                    bs_limo_display_resultsV4(pred,flags(flagIdx).group.folder,[],'MCC',5,'alpha',alpha,'style','topoPlotOnly')
                    be_print('file',[thisFileName '_topoOnly'],'prepare','plosBio','res','-r300','column',2)
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
%             continue
            try
                load(flags(flagIdx).subj(1).LIMO)
                tmpPred  = length(LIMO.design.nb_conditions)+pred*4;
                thisFileName = sprintf(filePath ,LIMO.design.XDesc{tmpPred},flags(flagIdx).name,tmpPred);
                
                
                bs_limo_display_resultsV4(100+pred,flags(flagIdx).group.folder,[],'MCC',5,'alpha',alpha,'interaction','yes')
                if g.print
                    be_print('file',thisFileName,'prepare','plosBio','res','-r150','column',2)
                    
                    bs_limo_display_resultsV4(100+pred,flags(flagIdx).group.folder,[],'MCC',5,'alpha',alpha,'interaction','yes','style','topoPlotOnly')
                    be_print('file',[thisFileName '_topoOnly'],'prepare','plosBio','res','-r300','column',2)
                end
                %             close all
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
    
end