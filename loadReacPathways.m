%%0.- READING DATA FROM EXCEL
%Import DG0f, chrV, activities (a_i)
xlsxFile = 'ButyrateOxidation.xlsx';
nSheet   = 'BuOxidPaths';
nRange = 'A1:DA100';

%% Load all the variable names and its values
[Values,  Strings]   = xlsread(xlsxFile, nSheet, nRange, 'basic');
diffColValStr = length(Strings(1,:)) - length(Values(1,:));

%% Identifiers of Blocks

%%  Path & Pathway
[rowPathList, colPathList]           = find(strcmp('Pathway List', Strings));
[rowPathSelection, colPathSelection] = find(strcmp('Selected Pathway', Strings));
[~, colFullSto] = find(strcmp('Full stoichiometry', Strings));

%Pathway Names
Pathway = Strings(rowPathList:rowPathSelection - 2, colPathList+1);
Path = Values(rowPathList:rowPathSelection - 2, (colPathList+2:colFullSto-1) - diffColValStr );   %Pathway Values      
PathwayNames = Strings(rowPathList:rowPathSelection - 2, colFullSto);

%% Mass Balance
[rowC, colC]          = find(strcmp('C', Strings));
[rowMet, colMet]      = find(strcmp('Met', Strings));
[lastRowSto, ~] = find(strcmp('Electron', Strings));

massBalValues = Values(rowC + 1 : lastRowSto, (colC:colMet) - diffColValStr);
massBalComp   = Strings(rowC, colC:colMet);


%% Thermodynamics blocks
% [rowDG, colDG]          = find(strcmp(['ΔGºf' char(10) '(kJ/mol)'], Strings));
% [rowHG, colHG]          = find(strcmp(['ΔHºf' char(10) '(kJ/mol)'], Strings));
% [rowS0, colS0]          = find(strcmp(['ΔSºf' char(10) '(kJ/mol K)'], Strings));
colDG = colC - 4;
colS0 = colC - 1;

thermoValues = Values(rowC + 1 : lastRowSto, (colDG:colS0) - diffColValStr);
thermoNames  = Strings(rowC, colDG:colS0);

%Fix ThermoNames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thermoNames = replace(thermoNames, 'Δ', 'D');
thermoNames = replace(thermoNames, 'º', '0');
thermoNames = regexprep(thermoNames,'\n+','');
thermoNames = replace(thermoNames, '(kJ/mol)', '');
thermoNames = replace(thermoNames, '(kJ/mol K)', '');
thermoNames = replace(thermoNames, '(J/mol K)', '');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stoichiometry
[rowLabels, colStoStart]   = find(strcmp('eD ID', Strings));
[rowProtTransloc,   ~] = find(strcmp('# Protons transloc', Strings));
[rowATP,   ~]    = find(strcmp('# ATP', Strings));
[~, colFullSto] = find(strcmp('Full stoichiometry', Strings));

reacM        = Values(rowLabels+1 : rowProtTransloc-1, (colStoStart +1:colFullSto-1)- diffColValStr);
reacHeadings = Strings(rowLabels, colStoStart +1: colFullSto-1);
reacLabels   = Strings(rowLabels - 1, colStoStart +1: colFullSto-1);


protTranslocM = Values(rowProtTransloc:rowProtTransloc+1, (colStoStart +1:colFullSto-1)- diffColValStr); 
n_ATPV        = Values(rowATP, (colStoStart +1:colFullSto-1)- diffColValStr);

%Components Ka
[~, col_Ka1]  = find(strcmp('pKa�1', Strings));
[~, col_Ka2]  = find(strcmp('pKa�2', Strings));

%Metabolite names, concentrations, charge, eD identification and pKa
[~, col_StM]  = find(strncmp(Strings,'ai / P',2));
[~, col_chrV] = find(strcmp(Strings,'Charge'),1);
[~, col_StNames] = find(strcmp(Strings,'Abbrev. Name'),1);
rngXValues = rowLabels+1:rowProtTransloc-1;

StNames = Strings(rngXValues, col_StNames);  
StM     = Values(rngXValues, col_StM - diffColValStr);
chrV    = Values(rngXValues, col_chrV - diffColValStr);
eD      = Values(rngXValues, colStoStart - diffColValStr);    %Identify electron donor of the reaction
pKa1V    = Values(rngXValues, col_Ka1 - diffColValStr); 
pKa2V    = Values(rngXValues, col_Ka2 - diffColValStr); 
% Substitute NaN values by 0
ka1V(isnan(pKa1V)) = 0;
ka2V(isnan(pKa2V)) = 0;

%Store in the structure of states
St.eD = eD;    % Store electron donor of the reaction defined in the spreadsheet
St.pKa1V = pKa1V;   St.pKa2V = pKa2V;

%Identifies the phase of the components (IN or OUT of the cell) and the permeability of them to cross the cell membrane
[~, col_compPhase]        = find(strcmp(Strings, 'Phase'));
[~, col_compPermeability] = find(strcmp(Strings, 'Permeable'));
[~, col_StType] = find(strcmp(Strings, 'CompLabel'));
[~, col_StPlot] = find(strcmp(Strings, 'StPlot'));

%**************************************************************************
%% DATA TREATMENT: Manipulates matrices values and string names to be used afterwards

%% 1.- Mass Balance Components **********************************************
% massBalComp{end} = 'S-CoA';                       %Adds the CoA as last mass balance component to C, H, N, P, S, etc.
massBalComp = strrep(massBalComp, '-', '_');      %Replaces the - by _ for the deprotonated components (For Ad- in this case, to avoid reading problems using the minus.)

%Assigns to each variable (C, H, N, etc.) the number of each component per each one of the states defined
%in St.StNames (e.g: C=3 for Propionate)
for i = 1:length(massBalComp)
    Reac.MassBal.(char(massBalComp(i))) = massBalValues(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.- STRINGS HANDLING  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Separates the headings from the rest of the matrix
% [rowReacLabel, ~]          = find(strcmp('ReacLabel', Strings)) ;   %Identify row where all the states start
% [rowHeadings, col_StNames] = find(strcmp('Abbrev. Name', Strings)) ;   %Identify row where all the states start
% [row_ATPV, ~]              = find(strcmp('# ATP', Strings));
% [row_protTranslocV, ~]     = find(strcmp('# Protons transloc', Strings));

% reacLabels   = Strings(rowReacLabel, :);  %Reaction Labels for electron Carriers (eC), fixed products (Fp) or Calculated products (Cp)
% reacHeadings = Strings(rowHeadings,  :);

%Collects the name that are NOT Headings to collect later on each group of names
% names = Strings(rowHeadings + 1 : row_protTranslocV - 1,:);

%Recalls the name of the variables
% StNames = names(1:end, col_StNames);

% %Identifies the phase of the components (IN or OUT of the cell) and the permeability of them to cross the cell membrane
% col_compPhase = strcmp(Strings(2,:),'Phase');
% col_compPermeability = strcmp(Strings(2,:),'Permeable');
% col_StType = strcmp(Strings(2,:),'CompLabel');
% col_StPlot = strcmp(Strings(2,:), 'StPlot');

St.compPhase  = strcmp(Strings(rngXValues, col_compPhase), 'IN');           %1 if it is IN and 0 if it is OUT
St.compPermeability = strcmp(Strings(rngXValues,col_compPermeability),  'Y');     %1 if it is PERMEABLE and 0 if it is NOT PERMEABLE
St.StType = (Strings(rngXValues, col_StType));                         %'m' for metabolite, 'eC' for electron Carrier, 'l' for loop component and 'r' for ratio-calculated component
St.StPlot = (Strings(rngXValues, col_StPlot));      %Select which variable to plot
% Identify positions of the variables to close the stoichiometry later
StNamesId = strrep(StNames,   '-', '_');      %Replaces the - by _ for the deprotonated components
StNamesId = strrep(StNamesId, '+',  '');     %Removes the + to avoid errors in the structure
StNamesId = strrep(StNamesId, ' ', '_');    %Replaces spaces for variables with _ (in case there is any)

%Differentiate between components IN and OUT of the cell
for i= 1:length(St.compPhase)
    if St.compPhase(i) == 0
        StNamesId(i) = strrep(StNames(i), StNames(i), strcat(StNamesId(i), 'out'));
    end
end

%Update St.StNames to differentiate IN and OUT Phases
St.StNames  = StNames;
St.StNames(St.compPhase == 0) = StNamesId(St.compPhase ==0);

%Assigns position identifier to each variable Name
for i = 1:length(StNames), St.id.(char(StNamesId(i))) = i; end

%**************************************************************************
%% 3.- NUMBERS RE-SCHEDULING

%Reschedule Values, Proton Translocation Matrix and number of ATP in a reaction
% diffRowStrVal = length(Strings(:,1)) - length(Values(:,1));
% diffColStrVal = length(Strings(2,:)) - length(Values(1,:));
% protTranslocM = Values(row_protTranslocV - diffRowStrVal : row_protTranslocV + 1 - diffRowStrVal, :);
% n_ATPV        = Values(row_ATPV - diffRowStrVal,:);
% recalc_DG     = Values(row_ATPV - diffRowStrVal +1,:);
% Values        = Values(1:row_protTranslocV - diffRowStrVal-1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Potential differential between inside and outside the cell membrane ### USER-EDITED ###
%(Could be read also from Excel)
Param.Rth =   0.00831;       Param.F   =   96.485;
Param.Max_Conc = 1e-2;       Param.Min_Conc = 1e-6;
Param.T   = 298.15;

% Adding Constants and Concentrations of the components
%Collect variables for DG0f, chrV and a_i
St.StM      = StM;           
Param.chrV  = chrV;
% Param.DG0f  = Values(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define DGf of the components Concentrations
Param.G0fM = thermoValues(:, 1);
Param.H0fM = thermoValues(:, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Energetics parameters subject to change
% Param.DG_ATP    = -60; % kJ/mol (hydrolysis of ATP)
% n_H_ATPsynth    =   9;
minProtTransloc =  -2;
maxProtTransloc =   6;
num_H_ATP =   3;

%Define range of proton translocation and its names for collecting results
%as per their net proton translocation
Param.netHTranslocRange = (minProtTransloc :(1/3): maxProtTransloc);
Param.netHTranslocRange = round(fliplr(Param.netHTranslocRange),2);

for i = 1:length(Param.netHTranslocRange)
    Param.netHTranslocNames{i} = sprintf('NetHTransloc_%s',num2str(Param.netHTranslocRange(i)));
end

Param.netHTranslocNames= strrep(Param.netHTranslocNames, '.', '_')';


%Change minus sign for extra underscore
Param.netHTranslocNames = strrep(Param.netHTranslocNames, '-', '_');


% %**************************************************************************
%% 4.- Loading Stoichiometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Finds the initial point of the reaction (Propionate transported in) and
%the end of the pathaway (Full Stoichiometry) (Note: Deducts the difference
%between values and strings to match columns between them)
Reac.id.Start  = colStoStart + 1 - diffColValStr;  % +1 because the starting of the full reaction is always 1 column after eD ID
Reac.id.End    = colFullSto  - 1 - diffColValStr ;

%Collects the whole potential stoichiometry for all possible pathways
%except the final full stoichiometry (-1)
% reacM = Values(:, Reac.id.Start:Reac.id.End-1);
Reac.FullStoM = zeros(length(reacM(:,1)),length(Pathway));

%Full stoichiometry from Excel (USED FOR LATER COMPARISONS ONLY)
% Reac.stoM_xlsx = Values(:, Reac.id.End);
Reac.stoM_xlsx =  Values(rowLabels+1 : rowProtTransloc-1,colFullSto - diffColValStr);

%Loads Reaction Pathway(Names) and Paths (Matrix with reactions orders)
Reac.Pathway = Pathway;         Reac.Path = Path;       Reac.PathwayNames = PathwayNames;
%Collects all the headings of the reactions specified in Excel
Reac.NamesV  = reacHeadings;
Reac.LabelsV = reacLabels;

%Converts to 1s and 0s those reactions that have been labeled in the Excel
%as membrane proteins
% Reac.MemProt = strcmp(Strings(end,(Reac.id.Start+DiffStrVal:Reac.id.End+DiffStrVal)),'Y');
Reac.protTranslocM = protTranslocM;
Reac.n_ATPV        = n_ATPV;

%**************************************************************************
numPathways = length(Pathway);
for i = 1:numPathways
    
    reacPath = Reac.Pathway(i);
    
    %Trim reactions, orders and names
    PathOrder = Path(i,:);
    Reac.(char(reacPath)).PathOrder    = PathOrder           (Path(i,:)   > 0);
    Reac.(char(reacPath)).Names        = Reac.NamesV         (Path(i,:)   > 0);
    Reac.(char(reacPath)).Labels       = Reac.LabelsV        (Path(i,:)   > 0);
    Reac.(char(reacPath)).protTransloc = Reac.protTranslocM  (:,Path(i,:) > 0);
    Reac.(char(reacPath)).n_ATP        = Reac.n_ATPV         (:,Path(i,:) > 0);
    Reac.(char(reacPath)).stoM_trimmed = reacM               (:,Path(i,:) > 0);
    
    %4b. Electron Carriers re-schedule
    %%
    %Definition of reduced e-Carriers (Could be read also from Excel by
    %labelling eC's)
    
    %Identify electron Carriers
    St_eC = St.StNames(strcmp(St.StType,'eC'));
    St_eB = St.StNames(strcmp(St.StType,'eB'));
    %     eC_finder = zeros(length(St_eC),1);
    eB_finder = zeros(length(St_eB),1);
    %     %Find reactions with electron carrier regeneration
    %     for z = 1:length(St_eC)
    %         eC_finder(z) = sum(strncmp(Reac.NamesV, St_eC(z)));
    %     end
    Param.eCarriers = St_eC;
    Param.eC_Cons   = St.StNames(strcmp(St.StType,'eCox'));
    Param.eC_Cons   = strrep(Param.eC_Cons, '+', '');
    
    
    %     Param.eC_Cons   = St_eC(eC_finder == 0);
    %     Param.eCarriers = {'F420red'};
    
    %     Param.eCarriers = {'NADH', 'NADPH', 'UQred', 'Fdred', 'FADH2', 'F420ox'};
    
    %Find reactions with electron bifurcation
    for z = 1:length(St_eB)
        eB_finder(z) = sum(strncmp(Reac.NamesV, St_eB(z),6));
    end
    
    Param.eBifurc = St_eB(eB_finder == 1);
    
    
    %     Param.eCarriers = {'CoM_S_S_CoB', 'F420red'};
    Param.Reox_reac = strcat(Param.eCarriers, '_reox');
    Param.eB_reac   = strcat(Param.eBifurc, '_eB');
    Param.eB_reac = strrep(Param.eB_reac, '-', '_');
    
    
    % e-Carriers reoxidations (or regeneration / moiety conservation)
    for j = 1:length(Param.eCarriers)
        %Identify reoxidation reactions
        %         Reac.id.(char(Param.Reox_reac(j)))  = find(strcmp(Reac.NamesV, Reac.NamesV(strncmp(Reac.NamesV,Param.eCarriers(j),4))));
        eC_ConsReax = strncmp(Reac.(char(Pathway(i))).Names,Param.eCarriers(j),4);
        
        %         if isempty(Reac.Names.(char(Pathway(i)))(strncmp(Reac.Names.(char(Pathway(i))),Param.eCarriers(j),4))),
        if sum(eC_ConsReax) == 0
            Reac.(char(Pathway(i))).id.(char(Param.Reox_reac(j))) = -1;
        else
            id_eC_Cons_Reax = find(eC_ConsReax == 1);
            Reac.(char(Pathway(i))).id.(char(Param.Reox_reac(j))) = id_eC_Cons_Reax;
            %             Reac.id.(char(Pathway(i))).(char(Param.Reox_reac(j))) = find(strcmp(Reac.Names.(char(Pathway(i))), Reac.Names.(char(Pathway(i)))(strncmp(Reac.Names.(char(Pathway(i))),Param.eCarriers(j),4))));
        end
    end
    
    id_eB_Reax = strcmp(Reac.LabelsV, 'eB');
    
    if any(id_eB_Reax == 1)
        for j = 1:length(Param.eBifurc)
            eB_Reax = strncmp(Reac.(char(Pathway(i))).Names, Param.eBifurc(j),5);
            
            id_eB_Reax = find(eB_Reax == 1);
            Reac.(char(reacPath)).id.(char(Param.eB_reac(j))) = id_eB_Reax;
        end
    end
    
    
    
    %%    %4c. Calculates the full stoichiometry, multiplying by 0 those columns that do
    
    %Trim stoichiometry and Sorts every reaction by column and puts them in the appropriate order
    %     Reac.(char(reacPath)).stoM_trimmed = Reac.(char(reacPath)).stoM(:,Path(i,:)>0);
    
    [col_sorted, col_order]                = sort(Reac.(char(reacPath)).PathOrder);
    Reac.(char(reacPath)).stoM_ord         = Reac.(char(reacPath)).stoM_trimmed(:,col_order);
    %    not belong to the Pathway highlighted
    Reac.(char(reacPath)).stoM  = bsxfun(@times, reacM, Path(i,:)>0);
    
    if numPathways > 1
        
        %Recalculates e-Carriers, H+ in and ATP consumption stoichiometries for each one of the defined Pathway.
        [recalcStoM_ord] = recalcStoich(Reac, St, reacPath, Param);
        
        Reac.(char(reacPath)).stoM_ord = recalcStoM_ord;
        
    end
    % 4b. Ordering Stoichiometries and eliminating columns that are 0
    % Load every stoichiometry IN ORDER and only with the reactions that are running
    
    %Trim stoichiometry and Sorts every reaction by column and puts them in the appropriate order
    %     Reac.(char(reacPath)).stoM_trimmed = Reac.(char(reacPath)).stoM(:,Path(i,:)>0);
    
    %     [col_sorted, col_order]                = sort(Reac.(char(reacPath)).PathOrder);
    %     Reac.(char(reacPath)).stoM_ord         = Reac.(char(reacPath)).stoM_trimmed(:,col_order);
    Reac.(char(reacPath)).Names_ord        = Reac.(char(reacPath)).Names(:,col_order);               %Write names of steps for each reaction
    Reac.(char(reacPath)).Labels_ord       = Reac.(char(reacPath)).Labels(:,col_order);               %Write names of Labels for each reaction in the correct order
    Reac.(char(reacPath)).protTransloc_ord = Reac.(char(reacPath)).protTransloc(:,col_order);
    Reac.(char(reacPath)).n_ATP_ord        = Reac.(char(reacPath)).n_ATP(:,col_order);
    
    %4d. Comprobation of Full Stoichiometry being the same for all the pathways
    Reac.FullStoM(:,i) = sum(Reac.(char(reacPath)).stoM_ord, 2);
    
    %4e. IDENTIFY ATP AT SLP
    % Num of ATP obtained in the pathway by SLP (do only for first pathway
    % since all reactions should yield same final stoichiometry
    if i == 1
        Param.n_ATP_SLP = sum(Reac.(char(reacPath)).n_ATP_ord); %Num of ATP produced per SLP
        n_H_ATP_SLP = num_H_ATP * Param.n_ATP_SLP;
    end
    
    %4e. Create combinations for proton translocations for each mode
    [allProtTranslocComb] = calcProtTranslocComb(Reac.(char(reacPath)).protTransloc_ord, minProtTransloc, maxProtTransloc, n_H_ATP_SLP);
    Reac.(char(reacPath)).AllProtTranslocComb = allProtTranslocComb;
    
    
    %Sum of total proton translocations (Watch out if it is the +3 or another value: MOST LIKELY DOES NOT NEED TO BE HERE)
    %     Reac.sumNetProtTranslocComb.(char(reacPath)) = sum(Reac.AllProtTranslocComb.(char(reacPath)),2) + 3;
    
    
end

% Check that all stoichiometries are the same
stoM_check = bsxfun(@minus, Reac.stoM_xlsx, Reac.FullStoM);

if all(stoM_check(:) == 0)
    sprintf('Every pathway yields the same final stoichiometry')
else
    j = find(all(stoM_check == 0)==0);
    sprintf('Reaction Pathway %d does not yield the same full stoichiometry calculated in Excel\n', j)
end

%% 4.- IDENTIFY SUBSTRATES AND PRODUCTS
for i = 1:length(Reac.Pathway)
    
    reacPath = Reac.Pathway(i);
    [constConcM, varConcM, posConcV, idCalcReax, numReac, solvePathwayID] = idProdsSubs(Reac, St, reacPath);
    Reac.(char(reacPath)).constConcM      = constConcM;
    Reac.(char(reacPath)).varConcM        = varConcM;
    Reac.(char(reacPath)).posConcV        = posConcV;
    Reac.(char(reacPath)).idCalcReax      = idCalcReax;
    Reac.(char(reacPath)).numReac         = numReac;
    Reac.(char(reacPath)).PathSolvingMode = solvePathwayID;
    
end

%% 5.- IDENTIFY LOOPS AND/OR ELECTRON BIFURCATION
[idLoop] = idLoops(St, Reac);


%% 6.- Checking Mass-Balances  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mass balance check for all stoichiometries
for i = 1:length(massBalComp)
    for j=1:length(Pathway)
        Reac.MassBal.Check(i,j) = sum(sum(bsxfun(@times, Reac.(char(Pathway(j))).stoM, Reac.MassBal.(char(massBalComp(i))))));
    end
end

%Mass Balance check of all reactions
if any(or(Reac.MassBal.Check > 1e-3, Reac.MassBal.Check < -1e-3)) == 0
    sprintf('Mass Balance is Correctly Checked')
else
    sprintf('Mass Balance is NOT PROPERLY CLOSED. PLEASE CHECK')
    %     sprintf('Mass Balance is NOT PROPERLY CLOSED FOR %s', char(MassBal_Comp((Reac.MassBal.Check ~= 0)))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 7.- Load Sensitivity Analysis Variables

%% Reading the variables from the sheet from the excel
nSheet   = 'SensAnalysis';
nRange = 'A1:DA100';

%% Load all the variable names and its values
[ValuesSens,  StringsSens]   = xlsread(xlsxFile, nSheet, nRange, 'basic');

%Find positions of combinations and of reference values
idrefVal       = find(strcmp(StringsSens, 'refVal'));
idvarNames     = find(strcmp(StringsSens, 'varNames'));
idvarUnits     = find(strcmp(StringsSens, 'Units'));
idSensAnalysis = find(strcmp(StringsSens, 'evalCombination'));

diffRowsStrVals = idrefVal - idSensAnalysis;

%Find variables used for sensitivity analysis
combNamesFull = StringsSens(idvarNames,:);  %Need still to trim
varUnitsFull  = StringsSens(idvarUnits,:);  

emptyChar = find(combNamesFull == "",1); %Identifies empty char
combNames = combNamesFull(2:emptyChar-1);
varUnits  = varUnitsFull(2:emptyChar-1); %Units of the variables

evalCombination = StringsSens(idSensAnalysis, 2);

%Find Reference Values for the sensitivity analysis
refVals = ValuesSens(idrefVal - diffRowsStrVals,:);
refVals(isnan(refVals)) = [];  %Trim reference values selected in Excel

%Find Values for each variable of the sensitivity analysis
FullValues = ValuesSens(idvarUnits-diffRowsStrVals+1:length(ValuesSens(:,1)), 1:length(refVals));

for i = 1:length(combNames)
    Combination.(char(combNames(i))) = FullValues(:,i);
    Combination.(char(combNames(i)))(isnan(Combination.(char(combNames(i))))) = [];
end

%Create reference combination vector
refComb = 0 * refVals;

for i = 1:length(refComb)
    refComb(i) = Combination.(char(combNames(i)))(refVals(i));
end

%Collect all data under Parameters structure
Param.refComb = refComb;
Param.sensAnalysis   = evalCombination;
Param.combNames      = combNames;
Param.combNamesUnits = varUnits;


%% Clears workspace
clearvars -except Reac Param St idLoop Combination

