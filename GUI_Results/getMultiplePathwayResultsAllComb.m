% clearvars -except PrintResults St Reac Results idLoop Param

% Need to explain CLEARLY still what is going on in this script. Consider
% also split into two scripts # MPG: 16/9/19

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  DATA COLLECTION/GENERATION   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear  and load all data
clear all
% load ResultsPathwayEvaluation.mat
load SimResults200scn.mat

%% 1. DEFINE PATHWAY TO PLOT THE RESULTS OBTAINED
refPoint = Param.refComb;
nPrevReac = 0;    %In case this plot is combined to another one

%%Obtain the pathway for the highest ATP Yield (Alternative, define a pathway to plot)
numPathways = length(Reac.Pathway);
maxNetATP   = zeros(numPathways, 1);


%Collect highest ATP Yield for EACH pathway
for i = 1:numPathways
    reacPath = Reac.Pathway(i);
    if isfield(PrintResults.FeasComb, reacPath)
        maxNetATP(i) = max(PrintResults.FeasComb.(char(reacPath)).netATPV);
    else
        maxNetATP(i) = -10;
    end
end

%Identify the pathway with the maximum ATP yield to plot it
% idMaxATP       = find(maxNetATP == max(maxNetATP));
% reacPathMaxATP = Reac.Pathway(idMaxATP(1));   %First pathway with largest ATP yield
% reacPathPlot   = reacPathMaxATP;                %Pathway to plot

%% 2. FILTER BY REFERENCE VALUE
for i=1:length(Reac.Pathway)
    
    
    reacPathPlot = Reac.Pathway(i);
    
    if isfield(PrintResults.FeasComb, reacPathPlot)
        %%Select a variable to sort results from
        % varToSort = 'H2';
        %Labels and Inputs results to plot Table ,
        InputParNames = St.combinationListNames;
        varToSort     = InputParNames(1);
        colSort       = find(strcmp(St.combinationListNames, varToSort));
        
        %Prepare text for legend
        outputVarNames  = {'\DeltaG_{ATP}', 'Ratio H^+ / ATP', '[CoA-SH]', 'pH_{in}', 'P_{H2}', '[CO_2]_{Diss}',   'T', 'pH_{out}'};
        unitsVarLabel   = {'kJ/mol'       ,                '',       'mM',        '',     'Pa',            'mM', '^oC',         ''};
        outputVarLegend = outputVarNames(colSort);
%         unitsVar        = unitsVarLabel(colSort);
         
        %Collect variables to be sorted by an input or net ATP
        Inputs           = PrintResults.FeasComb.(char(reacPathPlot)).Inputs;
        netATPV          = PrintResults.FeasComb.(char(reacPathPlot)).netATPV;
        ConcM            = PrintResults.FeasComb.(char(reacPathPlot)).ConcV;
        protTranslocComb = PrintResults.FeasComb.(char(reacPathPlot)).ProtTranslocComb;
        eC_Ratios_Res    = PrintResults.FeasComb.(char(reacPathPlot)).eC_Ratios_Res;
        DGrM_Reac        = PrintResults.FeasComb.(char(reacPathPlot)).DGrV;
        
        %Set the reference vector
        refValsV = Param.refComb(:,[1:colSort-1, colSort+1:end]);
        %Find all the lines from the feasible reactions that have the same equal parameters except the selected variable to be changed
        idRefVals = find(all(refValsV == Inputs(:,[1:colSort-1, colSort+1:end]),2));
        %Filter the results by the reference variables for inputs and net ATP
        Inputs_RefVals           = Inputs(idRefVals,:);
        netATP_RefVals           = netATPV(idRefVals,:);
        ConcM_RefVals            = ConcM(idRefVals,:);
        eC_Ratios_RefVals        = eC_Ratios_Res(idRefVals,:);
        protTranslocComb_RefVals = protTranslocComb(idRefVals,:);
        DGrV_RefVals             = DGrM_Reac(idRefVals,:);
        
        %Collect lengths of each of the variables
        length_Inputs           = length(refPoint);
        length_ConcM            = length(St.StM);
        length_protTranslocComb = length(protTranslocComb_RefVals(1,:));
        length_eC               = length(eC_Ratios_RefVals(1,:));
        length_DGrV             = length(DGrV_RefVals(1,:));
        
        %Combine all the results to be sorted
        ResultsToSort = [Inputs_RefVals, ConcM_RefVals, eC_Ratios_RefVals, protTranslocComb_RefVals, DGrV_RefVals, netATP_RefVals];
        
        %% Sort Results
        sorting = 'refVals';
        
        if strcmp(sorting, 'refVals')
            sortedResults      = sortrows(ResultsToSort, colSort, 'ascend');
            %Trim sortedResults into their own variables
            InputsSorted       = sortedResults(:,1:length_Inputs);
            ConcM_Sorted       = sortedResults(:,length_Inputs + 1:length_Inputs + length_ConcM);
            eC_Ratios_Sorted   = sortedResults(:,length_Inputs + length_ConcM + 1:length_Inputs+length_ConcM+length_eC);  %TO CORRECT STILL
            protTranslocSorted = sortedResults(:,length_Inputs + length_ConcM + length_eC + 1:length_Inputs + length_ConcM + length_eC + length_protTranslocComb);
            DGrM_Sorted        = sortedResults(:,length_Inputs + length_ConcM + length_eC + length_protTranslocComb + 1:end-1);
            netATPSorted       = sortedResults(:,end);
        end
        
%         %Collecting Results in Structure
%         plotResults.(char(reacPathPlot)).InputsSorted       = InputsSorted;
%         plotResults.(char(reacPathPlot)).ConcM_Sorted       = ConcM_Sorted;
%         plotResults.(char(reacPathPlot)).eC_Ratios_Sorted   = eC_Ratios_Sorted;
%         plotResults.(char(reacPathPlot)).protTranslocSorted = protTranslocSorted;
%         plotResults.(char(reacPathPlot)).DGrV_Sorted        = DGrV_Sorted;
%         plotResults.(char(reacPathPlot)).netATPSorted       = netATPSorted;
        
%         %Identify each singular value used to sort each variable
%         valsSorted = unique(InputsSorted(:,colSort));
%         % valsH2 = sort(valsH2, 'descend');
%         idValsSorted = zeros(length(valsSorted),1);
%         for k = 1:length(valsSorted)
%             idValsSorted(k) = find(InputsSorted(:,colSort) == valsSorted(k), 1);
%         end
        
        %         %% Parameters for units conversion units
        %         Cels2Kelv = 273.15;    DG_H2 = 18.3;
        
        %% Define the number of results to be plotted
        %         idPosValsSorted  = idValsSorted(netATPSorted(idValsSorted) >0);  %Plot only positive net ATPs
        %         numRes = length(idPosValsSorted);
        numLoops  = idLoop.(char(reacPathPlot)).numLoops;
                
        %% PATHWAY DATA
        plotSt = strcmp(St.StPlot, 'p');  %Only the metabolites labelled as 'p' are plotted       

        %Number of reactions that are not eB nor eC reoxidations
        numRc_eB = sum(strcmp(Reac.(char(reacPathPlot)).Labels_ord, 'eB'));
        numRc_eC = sum(strcmp(Reac.(char(reacPathPlot)).Labels_ord, 'eC'));
        RcLength = length(Reac.(char(reacPathPlot)).Names_ord);
        nPlotRcNames = RcLength - numRc_eB - numRc_eC;
        
        %Collect  Names for XTicks
        reacNames =  Reac.(char(reacPathPlot)).Names_ord(1:nPlotRcNames);        
        barLineWidth = 1;
        
        %Remove all metabolites that do not participate in any reaction
        stoM   = Reac.(char(reacPathPlot)).stoM_ord;    % Load stoichiometry of the pathway
        rows0  = find(all(stoM==0,2));                  % Identify the rows with metabolites not participating 
        plotSt(rows0) = [];                             % Delete all the rows that are zero
        StNames   = St.StNames;                         % Load the names of the metabolites
        
        
        %Arrange the values for Proton Translocations and ATP via Substrate
        %Level Phosphorylation
        protTranslocReac = protTranslocSorted(:,1:length(reacNames));                   %Proton translocation for the pathway reactions
        protTransloc_eC  = protTranslocSorted(:,length(reacNames)+1:end);               %Proton translocation for the electron carrier reoxidations (not ordered yet with the names)
        ATP_SLP_V           = Reac.(char(reacPathPlot)).n_ATP_ord(1:length(reacNames));    %ATP at substrate level phosphorylation reactions
                
        %Arrange the values for Gibbs Free Energies
        DGrM_Reac          = DGrM_Sorted(:,1:length(reacNames));
        DGrM_eC       = DGrM_Sorted(:,length(reacNames)+1:end);                         %Gibbs energy for the electron carrier reoxidations (not ordered yet with the names)        
        r_H_ATP_pos   = strcmp(InputParNames, 'ratio_H_ATP');
        ratio_H_ATP   = refPoint(1,r_H_ATP_pos);
        
        % Write Abbreviated names (R_i) for Reac Labels plot
        reacAbbrevNames = cell(size(Reac.NamesV));
        for l=1:length(Reac.NamesV)
            reacAbbrevNames{l} = sprintf('R_{%d}',l+nPrevReac);
        end
        %Identify which reactions from Pathway P_i
        reacID = zeros(1,length(reacNames));
        for j = 1:length(reacNames)
            reacID(j) = find(strcmp(reacNames(j), Reac.NamesV));
        end
        reacID = sort(reacID);
        reacPlotAbbrevNames = reacAbbrevNames(reacID);
        
        %In case of loop reactions, the loop reaction is labelled twice! 
        %Collect  Names for XTicks in case of loop reactions, proton
        %translocations, ATP at SL and calculate de Gibbs free energy
        if numLoops > 0
            for j = 1:numLoops
                loopName      = sprintf('Loop_%d', j);
                startLoop     = idLoop.(char(reacPathPlot)).(char(loopName)).Start;
                endLoop       = idLoop.(char(reacPathPlot)).(char(loopName)).End;
                reacNames     = [reacNames(1:endLoop), reacNames(startLoop), reacNames(endLoop+1:end)];
                reacPlotAbbrevNames   = [reacPlotAbbrevNames(1:endLoop), reacPlotAbbrevNames(startLoop), reacPlotAbbrevNames(endLoop+1:end)];
                DGrM_Reac     = [DGrM_Reac(:,1:endLoop),      DGrM_Reac(:,startLoop),      DGrM_Reac(:,endLoop+1:end)];
                ATP_SLP_V     = [ATP_SLP_V(:,1:endLoop),    ATP_SLP_V(:,startLoop),    ATP_SLP_V(:,endLoop+1:end)];
                protTranslocReac = [protTranslocReac(:,1:endLoop), protTranslocReac(:,startLoop), protTranslocReac(:,endLoop+1:end)];
            end
        end
                
        %Prepare metabolites concentrations to be plotted
        ConcM_Sorted(:,rows0) = [];
        plotConcM   = ConcM_Sorted(:,plotSt);     
        %Prepare concentrations names of the metabolites to be plotted
        StNames(rows0) = [];
        StNamesTrimmed = StNames(plotSt);         
%         reacNum = 1:length(reacNames);
        
        %THIS IS A TEMPORARY FIX (Consider possibility of labeling initial
        %metabolite in the excel)
        id_Start_Plot = find(strcmp(StNamesTrimmed, 'Pro_out'));   %Find the position of the metabolite where the plot should start
        startNamePlot = StNamesTrimmed(id_Start_Plot);             %Find the metabolite name in the StNames vector
        StNamesTrimmed(id_Start_Plot) = [];                        %Remove the metabolite from their original location in the StNames vector
        plotStNames = [startNamePlot; StNamesTrimmed];             %Put the metabolite selected as the first name in the plotNames vector
        plotStNames = strrep(plotStNames, '_out', ' (out)');       %Replace the label _out by " (out)
        
        %Apply the same change to the concentration matrix than to the
        %names concentration matrix
        startConcPlot = plotConcM(:,id_Start_Plot);
        plotConcM(:,id_Start_Plot) = [];
        plotConcM = [startConcPlot, plotConcM]';
                          
        

        
        
        %% Prepare eC ratios to plot        
        if ~isempty(Param.eB_reac)
            eB_carrier    = 'CoB-SH';  %Select between 'CoB-SH' or 'Fdred'
        else
            eB_carrier = '';
        end
        
        %Calculate ratio of electron bifurcation  (THIS PART NEEDS STILL TO
        %BE FIXED)
        if ~isempty(eB_carrier)
            switch eB_carrier
                case 'CoB-SH'
                    eB_Ratios = ConcM_Sorted(St.id.CoB_SH,:) ./ ConcM_Sorted(St.id.CoM_S_S_CoB,:);
                    ratioNames = {'CoB-SH/CoM-S-S-CoB', 'F420_{ox}/F420_{red}'};
                case 'Fdred'
                    eB_Ratios = ConcM_Sorted(St.id.Fdred,:) / ConcM_Sorted(St.id.Fdox,:);
                    ratioNames = {'Fd_{red}/Fd_{ox}', 'F420_{ox}/F420_{red}'};
            end
            plotRatios = [eB_Ratios', eC_Ratios_Sorted];
            reoxeCNames = Results.Combination_1.(char(reacPathPlot)).eC_ReoxNames;   %any combination would work
            reoxeCNames = strrep(reoxeCNames, '_reox', '_{reox}');
            reacReoxAbbrevNames = {'e-Bifurc',char(reoxeCNames)};
            
            %Identify position of electron carrier and electron bifurcator
            id_eC_R = find(strcmp(Reac.(char(reacPathPlot)).Labels_ord, 'eC'));
            id_eB_R = find(strcmp(Reac.(char(reacPathPlot)).Labels_ord, 'eB'));
            
            %Number of protons translocated for electron carriers and
            %for reaction with electron bifurcation
            protTransloc_eC = protTranslocSorted(:, id_eC_R);
            protTransloc_eB = protTranslocSorted(:, id_eB_R);
            plotHTransloc   = [protTransloc_eB, protTransloc_eC];
            
            %Names of reactions for XTickLabels
            reoxReacNames = [Reac.(char(reacPathPlot)).Names_ord(id_eB_R), Reac.(char(reacPathPlot)).Names_ord(id_eC_R)];
            
        else
            ratioNames = Results.Combination_1.(char(reacPathPlot)).eC_ReoxNames;
            reacReoxAbbrevNames = strrep(ratioNames, '_reox', '_{reox}');
            ratioNames = strrep(ratioNames, '_reox', '');
            
            eCarrierRatioLabel = cell(length(ratioNames),1);
            eCarrierLabel = eCarrierRatioLabel;
            for z = 1:length(Param.eCarriers)
                
                idReox_Labels = find(strcmp(Param.eCarriers(z), ratioNames));
                if ~isempty(idReox_Labels)
                    %         reoxLabels(z) = strcat(Param.eCarriers(idReox_Labels),'/',Param.eC_Cons(idReox_Labels));
                    eCarrierLabel(idReox_Labels) = Param.eCarriers(z);
                    eCarrierRatioLabel(idReox_Labels) = strcat(Param.eCarriers(z),'/',Param.eC_Cons(z));
                end
            end
            plotRatios = eC_Ratios_Sorted;
        end
                
        %Prepare proton translocation reorderred and Gibbs free energies
        %according to eC reoxidation lables
        lastNonEcReac = find(strcmp(Reac.(char(reacPathPlot)).Labels, 'Fp')); %Lsat reaction before electron carriers
        
        eC_Reox_Names = Reac.(char(reacPathPlot)).Names_ord(lastNonEcReac+1:end);
        reoxReacID = zeros(1, length(eC_Reox_Names));
        for k = 1:length(eC_Reox_Names)
            reoxReacID(k) = find(strncmp(eC_Reox_Names(k), reacReoxAbbrevNames, 4));
        end
       
        %Reorder the columns accordingly to each electron Carrier
        plotProtTransloc_eC = protTransloc_eC(:,reoxReacID);
        plotDGrM_eC         = DGrM_eC(:,reoxReacID);
                  
        %Collecting Results in Structure
        plotResults.(char(reacPathPlot)).InputsSorted         = InputsSorted;      
        plotResults.(char(reacPathPlot)).netATPSorted         = netATPSorted;      

        plotResults.(char(reacPathPlot)).plotConcM            = plotConcM;
        plotResults.(char(reacPathPlot)).eC_Ratios_Sorted     = eC_Ratios_Sorted;
        plotResults.(char(reacPathPlot)).plotprotTranslocReac = protTranslocReac;
        plotResults.(char(reacPathPlot)).plotprotTransloc_eC  = plotProtTransloc_eC;
        plotResults.(char(reacPathPlot)).ATP_SLP_V            = ATP_SLP_V;
        plotResults.(char(reacPathPlot)).DGrM_Reac            = DGrM_Reac;
        plotResults.(char(reacPathPlot)).DGrM_eC              = plotDGrM_eC;
        
        %Collect labelings in structures
        plotResults.(char(reacPathPlot)).plotStNames          = plotStNames;
        plotResults.(char(reacPathPlot)).eCarrierRatioLabel   = eCarrierRatioLabel;        
        plotResults.(char(reacPathPlot)).reacPlotAbbrevNames  = reacPlotAbbrevNames;
        plotResults.(char(reacPathPlot)).reacFullNames        = reacNames;
        plotResults.(char(reacPathPlot)).reacReoxAbbrevNames  = reacReoxAbbrevNames;
        
    end
end

clearvars -except Reac Param idLoop St Results PrintResults Combination plotResults

