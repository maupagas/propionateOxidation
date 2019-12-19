% clearvars -except PrintResults St Reac Results idLoop Param

% Need to explain CLEARLY still what is going on in this script. Consider
% also split into two scripts # MPG: 16/9/19

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  DATA COLLECTION/GENERATION   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
idMaxATP       = find(maxNetATP == max(maxNetATP));
reacPathMaxATP = Reac.Pathway(idMaxATP(1));   %First pathway with largest ATP yield
reacPathPlot   = reacPathMaxATP;                %Pathway to plot

%% 2. FILTER BY REFERENCE VALUE

%%Select a variable to sort results from
% varToSort = 'H2';
%Labels and Inputs results to plot Table ,
InputParNames = St.combinationListNames;
varToSort     = InputParNames(5);
colSort       = find(strcmp(St.combinationListNames, varToSort));

%Prepare text for legend
outputVarNames  = {'\DeltaG_{ATP}', 'Ratio H^+ / ATP', '[CoA-SH]', 'pH_{in}', 'P_{H2}', '[CO_2]_{Diss}',   'T', 'pH_{out}'};
unitsVarLabel   = {'kJ/mol'       ,                '',       'mM',        '',     'Pa',            'mM', '^oC',         ''};
outputVarLegend = outputVarNames(colSort);
unitsVar        = unitsVarLabel(colSort); 

%Collect variables to be sorted by an input or net ATP
Inputs           = PrintResults.FeasComb.(char(reacPathPlot)).Inputs;
netATPV          = PrintResults.FeasComb.(char(reacPathPlot)).netATPV;
ConcM            = PrintResults.FeasComb.(char(reacPathPlot)).ConcV;       
protTranslocComb = PrintResults.FeasComb.(char(reacPathPlot)).ProtTranslocComb;
eC_Ratios_Res    = PrintResults.FeasComb.(char(reacPathPlot)).eC_Ratios_Res;
DGrV             = PrintResults.FeasComb.(char(reacPathPlot)).DGrV;

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
DGrV_RefVals             = DGrV(idRefVals,:);

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
        ConcM_Sorted       = sortedResults(:,length_Inputs+1:length_Inputs+length_ConcM);
        eC_Ratios_Sorted   = sortedResults(:,length_Inputs+length_ConcM+1:length_Inputs+length_ConcM+length_eC);  %TO CORRECT STILL
        protTranslocSorted = sortedResults(:,length_Inputs+length_ConcM+length_eC+1:length_Inputs+length_ConcM+length_eC+length_protTranslocComb);
        DGrV_Sorted        = sortedResults(:,length_Inputs+length_ConcM+length_eC+length_protTranslocComb+1:end-1);
        netATPSorted       = sortedResults(:,end);
end

%Identify each singular value used to sort each variable
valsSorted = unique(InputsSorted(:,colSort));
% valsH2 = sort(valsH2, 'descend');
idValsSorted = zeros(length(valsSorted),1);
for k = 1:length(valsSorted)
    idValsSorted(k) = find(InputsSorted(:,colSort) == valsSorted(k), 1);
end


%% Define the number of results to be plotted
idPosValsSorted  = idValsSorted(netATPSorted(idValsSorted) >0);  %Plot only positive net ATPs
numRes = length(idPosValsSorted);
numLoops  = idLoop.(char(reacPathPlot)).numLoops; 

%% Parameters for units conversion units
Cels2Kelv = 273.15;    DG_H2 = 18.3;

%% PATHWAY DATA
plotSt = strcmp(St.StPlot, 'p');  %Only the metabolites labelled as 'p' are plotted

%Remove all metabolites that do not participate in any reaction
stoM   = Reac.(char(reacPathPlot)).stoM_ord;
rows0  = find(all(stoM==0,2)); % all zeros
plotSt(rows0) = []; %Delete all the rows that are zero

% Collect concentration of metabolites of the selected results
ConcM_Res = ConcM_Sorted(idPosValsSorted,:);   
StNames   = St.StNames;

%Number of reactions that are not eB nor eC reoxidations
numRc_eB = sum(strcmp(Reac.(char(reacPathPlot)).Labels_ord, 'eB'));
numRc_eC = sum(strcmp(Reac.(char(reacPathPlot)).Labels_ord, 'eC'));
RcLength = length(Reac.(char(reacPathPlot)).Names_ord);
nPlotRcNames = RcLength - numRc_eB - numRc_eC;

%Collect  Names for XTicks
reacNames =  Reac.(char(reacPathPlot)).Names_ord(1:nPlotRcNames);
barLineWidth = 1;

%%% Find Abbreviated names for Reac Labels
reacAbbrevNames = cell(size(Reac.NamesV));
for i=1:length(Reac.NamesV)
    reacAbbrevNames{i} = sprintf('R_{%d}',i+nPrevReac);
end

reacID = zeros(1,length(reacNames));
for j = 1:length(reacNames)
    reacID(j) = find(strcmp(reacNames(j), Reac.NamesV));
end
reacID = sort(reacID);
reacPlotAbbrevNames = reacAbbrevNames(reacID);

%Assign translocated protons (ALL OF THEM NEED TO BE SORTED PROPERLY*****)
protTranslocV = protTranslocSorted(idPosValsSorted,:);
DGrV          = DGrV_Sorted(idPosValsSorted,1:length(reacNames));
DGrV_eC       = DGrV_Sorted(idPosValsSorted,length(reacNames)+1:end);
n_ATPV        = Reac.(char(reacPathPlot)).n_ATP_ord; 
r_H_ATP_pos   = strcmp(InputParNames, 'ratio_H_ATP');
ratio_H_ATP   = refPoint(1,r_H_ATP_pos);

protTranslocV = protTranslocV(:, 1:length(reacNames));
n_ATPV        = n_ATPV(1:length(reacNames));

if ~isempty(Param.eB_reac)
    eB_carrier    = 'CoB-SH';  %Select between 'CoB-SH' or 'Fdred'
else 
    eB_carrier = '';
end

%Collect  Names for XTicks
if numLoops > 0 
    for i = 1:numLoops
        loopName      = sprintf('Loop_%d', numLoops);
        startLoop     = idLoop.(char(reacPathPlot)).(char(loopName)).Start;
        endLoop       = idLoop.(char(reacPathPlot)).(char(loopName)).End;
        reacNames     = [reacNames(1:endLoop), reacNames(startLoop), reacNames(endLoop+1:end)];
        reacPlotAbbrevNames   = [reacPlotAbbrevNames(1:endLoop), reacPlotAbbrevNames(startLoop), reacPlotAbbrevNames(endLoop+1:end)];
        DGrV          = [DGrV(:,1:endLoop),      DGrV(:,startLoop),      DGrV(:,endLoop+1:end)];
        n_ATPV        = [n_ATPV(:,1:endLoop),    n_ATPV(:,startLoop),    n_ATPV(:,endLoop+1:end)];
        protTranslocV = [protTranslocV(:,1:endLoop), protTranslocV(:,startLoop), protTranslocV(:,endLoop+1:end)];
    end
end
posIndexes    = protTranslocV >= 0;


%Trim to """"plottable""" results
ConcM_Res(:,rows0) = [];
plotConcV   = ConcM_Res(:,plotSt);
StNames(rows0) = [];
StNames = StNames(plotSt);
reacNum = 1:length(reacNames);

%THIS IS A TEMPORARY FIX
id_Start_Plot = find(strcmp(StNames, 'Pro_out'));

startNamePlot = StNames(id_Start_Plot);
StNames(id_Start_Plot) = [];
plotStNames = [startNamePlot; StNames];
plotStNames = strrep(plotStNames, '_out', ' (out)');

startConcPlot = plotConcV(:,id_Start_Plot);
plotConcV(:,id_Start_Plot) = [];
plotConcV = [startConcPlot, plotConcV]';

%% eC REOXIDATION

%Get ratios of electron carriers 
eC_Ratios = eC_Ratios_Sorted(idPosValsSorted,:);

%Get concentrations result vector 
ConcM_Res = ConcM_Sorted(idPosValsSorted,:)';

%Calculate ratio of electron bifurcation
if ~isempty(eB_carrier)
    switch eB_carrier
        case 'CoB-SH'
            eB_Ratios = ConcM_Res(St.id.CoB_SH,:) ./ ConcM_Res(St.id.CoM_S_S_CoB,:);
            ratioNames = {'CoB-SH/CoM-S-S-CoB', 'F420_{ox}/F420_{red}'};
        case 'Fdred'
            eB_Ratios = ConcM_Res(St.id.Fdred,:) / ConcM_Res(St.id.Fdox,:);
            ratioNames = {'Fd_{red}/Fd_{ox}', 'F420_{ox}/F420_{red}'};
    end
    plotRatios = [eB_Ratios', eC_Ratios];
    reoxeCNames = Results.Combination_1.(char(reacPathPlot)).eC_ReoxNames;   %any combination would work 
    reoxeCNames = strrep(reoxeCNames, '_reox', '_{reox}');
    reacReoxAbbrevNames = {'e-Bifurc',char(reoxeCNames)};
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
    plotRatios = eC_Ratios;
end
   
%Identify position of electron carrier and electron bifurcator
id_eC_R = find(strcmp(Reac.(char(reacPathPlot)).Labels_ord, 'eC'));
id_eB_R = find(strcmp(Reac.(char(reacPathPlot)).Labels_ord, 'eB'));

%Number of protons translocated for electron carriers and 
%for reaction with electron bifurcation
protTransloc_eC = protTranslocSorted(idPosValsSorted, id_eC_R);
protTransloc_eB = protTranslocSorted(idPosValsSorted, id_eB_R);
plotHTransloc   = [protTransloc_eB, protTransloc_eC];

%Names of reactions for XTickLabels
reoxReacNames = [Reac.(char(reacPathPlot)).Names_ord(id_eB_R), Reac.(char(reacPathPlot)).Names_ord(id_eC_R)];


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% PLOT STARTS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PLOT DEFINITIONS
FontSize = 24;                  FontLegend = 20;
LabelRotation = 30;
deltaProt = 1;                      deltaDG = 5;  %# num of proton translocated per grid
minConc = log10(Param.Min_Conc);    maxConc = log10(Param.Max_Conc);
minProtTransloc = -3;               maxProtTransloc = 3;
nProt = minProtTransloc:maxProtTransloc;     %Range of protons translocated to plot
minRatio = -4;                      maxRatio = 4;
reacNameLabels = 0;   %0 for abbreviation (R1, R2, R3.... etc), 1 for fullname
reacLabelTicks = 0;   %1 to add ticks for reaction labels or 0 for not adding them

%Define colors of bars and markers
colorCode = 1;   % 0 for black and white, 1 for colored
markerList    = {'o',  's', '^',  'd', '*',  'x', 'p',  'h','<',  '>'};
lineStyleList = {'-', '--', ':', '-.', '-', '--', ':', '-.','-', '--'};

     
%Define bar colors 
if colorCode == 0
    posColor = [1  1  1];    %Color for positive values
    negColor = [1  1  1]; %Color for negative values
    ATPColor = [1  1  1];
    DGdissColor = [1 1 1];
    outerLimColor = [0.2  0.2  0.2 0.1];
    barTransparency = 0.7;
else
%     posColor = [0     0.5772  0];    %Color for positive values
%     negColor = [0.64  0.08    0.18]; %Color for negative values
    posColor      = [0.3922    0.8235    0.5294];    %Color for positive values
    negColor      = [1.0000    0.7059    0.7059]; %Color for negative values
    ATPColor      = [0.1176    0.4000    0.2118]; 
    outerLimColor = [0.64  0.08  0.18 0.1];
    DGdissColor   = [0.93 0.69 0.13];
    barTransparency = 0.8;
    barWidth = 1;
end

%Prepare strings to be printed in table and in legend
 colRatio = find(strcmp(InputParNames,'ratio_H_ATP'));
    ResData = [InputsSorted(idPosValsSorted,:), netATPSorted(idPosValsSorted,:)];
    
    %Modify strings in ResData (**************TO EVALUATE)
    numData = numel(ResData);
    strData = cell(size(ResData));
    %Re-structure dataLabels
    for i = 1:length(ResData(:,1))
        for j = 1:length(ResData(1,:))
            if abs(ResData(i,j)) >= 1
                if mod(ResData(i,j), 1) == 0 && j ~= colRatio
                    strData{i,j} = sprintf('%.0f', ResData(i,j));
                elseif j == colRatio
                    strData{i,j} = char(rationalTickLabels(ResData(i,j),3,0));
                else
                    strData{i,j} = sprintf('%.2f', ResData(i,j));
                end
            elseif ResData(i,j)<1 && ResData(i,j)>=0.1
                strData{i,j} = sprintf('%.2f', ResData(i,j));
            elseif ResData(i,j)<0.1 && ResData(i,j)>=0.001
                strData{i,j} = sprintf('%.0f', ResData(i,j) * 1000);
            elseif ResData(i,j)<0.001 && ResData(i,j)>=1e-6
                strData{i,j} = sprintf('%.3f', ResData(i,j) * 1000);
            elseif ResData(i,j)<1e-6 && ResData(i,j)>=0
                P_H2 = ResData(i,j) /(exp(-DG_H2/(Param.Rth * (Cels2Kelv +35)))) * 1e5;    %Pa
                strData{i,j} = sprintf('%.2f', P_H2);
            end
        end
    end


%Create legend labels
rnamesLegend = cell(1, length(idPosValsSorted));
rnamesTable  = rnamesLegend;
for i = 1:length(idPosValsSorted)
    rnamesTable{i} = sprintf('<html><tr align=center><col width=80><font size=+1>Config<sub>%d</sub>',i);
    rnamesLegend1  = sprintf('%s = %s %s', char(outputVarLegend), char(strData(i,colSort)), char(unitsVar));
    rnamesLegend2  = sprintf('Net ATP = %s', char(strData(i,end)));
    
    rnamesLegend{i} = sprintf('%s\n%s', char(rnamesLegend1), char(rnamesLegend2));
end

    

%Opens figure in Full Screen 
f = figure('units','normalized','outerposition',[0 0 1 1]);

%% PLOT DIMENSIONS: Main subplots in normalized units(Pathway and e-Carriers)
left_pos     = 0.10;
delta_vertical = 0.00;
bottom_posDG = 0.12; %Position of the DGPlot
width_1 = 0.65;      %Width of Pathway Plot
width_2 = 0.15;      %Width of e-Carriers Plot
height = 0.35;       %Height main subplots
height_HTransloc = 0.3;      %Height DG Plot
height_DG = 0.12;      %Height DG Plot
plot_space = 0.005;   %Separation between plots
bottom_protTransloc = bottom_posDG + height_DG + delta_vertical;
bottom_pos   = bottom_posDG + height_DG + height_HTransloc + delta_vertical*2; %Bottom position of the main plots

%% ************************************************
%************ TABLE DIMENSIONS ******************
%************************************************
isTable = 0;

if isTable == 1
     
    %%Write cnames as HTML code for table usage
    cnames2 = { '<html><font size=+1><p>&#916G<sub>ATP</sub> (kJ/mol)';         %DG_ATP
        %             '<html><font size=+1><tr align = "top">Ratio H<sup>+</sup>/ATP';   %ratio H/ATP
        '<html><font size=+1>Ratio H+/ATP';                    %ratio H/ATP
        '<html><font size=+1>CoA-SH (mM)</sub>';              %CoA-SH
        '<html><font size=+1>pH<sub>in</sub>';                 %pH in
        '<html><font size=+1>P<sub>H2</sub> (ppm)';            %Hydrogen Pressure
        '<html><font size=+1>CO<sub>2</sub> (mM)';              %CO2
        '<html><font size=+1>T (K)';                           %Temperature
        '<html><font size=+1>pH<sub>out</sub>';                %pH out
        '<html><font size=+1>Net ATP';};                       %ATP
    %Checkup for handmade table headers
    if length(cnames2) ~= length(St.combinationListNames)
        sprintf('It seems to be printing another pathway');
    end
    
    rnamesTable = cell(length(idPosValsSorted),1);

    % rnamesLegend = rnamesTable;
    
    % rnamesLegend1  = strcat('P_{H2} = ', strData(:,col_H2),'ppm');
    % rnamesLegend2 = strcat('Net ATP = ', strData(:,end));
    
    % rnamesLegend = sprintf('%s\n%s', char(rnamesLegend1), char(rnamesLegend2));
    
    % ResData =  PrintResults.FeasComb.(char(reacPath)).Inputs(numResults,:);
    % netATP  =  PrintResults.FeasComb.(char(reacPath)).netATPV(numResults,:);
   
    
    
    t_left   = 0.07;
    t_bottom = 0.82;
    t_width  = 0.095 * length(ResData(1,:));
    
    if length(idPosValsSorted) >= 3
        t_height = 0.06 * length(idPosValsSorted);
    elseif length(idPosValsSorted) == 2
        t_height = 0.07 * length(idPosValsSorted);
    elseif length(idPosValsSorted) == 1
        t_height = 0.1 * length(idPosValsSorted);
    end
    
    HeaderFontSize = 1;   HeaderWidth = 10;   columnWidth = 115;
    strCenteredData = strcat(sprintf('<html><tr align=center><td width=%d>', columnWidth), strData(:,:));
    InputParNames = strcat(sprintf('<html><font size=+%d>',HeaderFontSize), InputParNames(:));
    % rnamesHeader = strcat(sprintf('<html><font size=+%d>',HeaderFontSize), rnamesTable(:));
    rnamesHeader = rnamesTable;
    
    % Create the uitable to put data from INPUTS
    t = uitable(f,'Data',strCenteredData, 'ColumnName', cnames2, 'RowName',rnamesHeader,...
        'ColumnWidth',{columnWidth}, 'Units', 'normalized',...
        'Position', [t_left t_bottom t_width t_height]);
    t.FontSize = 14;
    
    %%% RESIZE THE FIRST COLUMN
    jscroll=findjobj(t);
    rowHeaderViewport=jscroll.getComponent(4);
    rowHeader=rowHeaderViewport.getComponent(0);
    % heightTable =rowHeader.getSize;
    rowHeader.setSize(80,360)
    %resize the row header
    newWidth = columnWidth;
    rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth, 0));
    heightTable = rowHeader.getHeight;
    rowHeader.setPreferredSize(java.awt.Dimension(newWidth, heightTable));
    rowHeader.setSize(newWidth, heightTable);
end

%% SUBPLOT 1: PLOT WHOLE PATHWAY
s1 = subplot(3,2,1);

%Plot concentrations for each metabolite in the pathway
concPlot = plot(log10(plotConcV));
ylabel('log_{10} C_{Metabolite}')

hL = legend(rnamesLegend, 'Orientation', 'horizontal', 'Location','SouthEast', 'box', 'off', 'FontSize', FontLegend);
hL.Position = [0.12 0.575 0.35 0.08];
%Set marker colours and sizes
for i = 1:numRes
    concPlot(i).MarkerSize = 10;
    concPlot(i).Marker = markerList{i};
    concPlot(i).LineStyle = lineStyleList{i};
    concPlot(i).MarkerFaceColor = [0.8 0.8 0.8];
    concPlot(i).MarkerEdgeColor = [0  0  0];
    concPlot(i).LineWidth = 1;
    concPlot(i).Color = [0 0 0];
end

%Sets characteristics of Top X-axis (Concentrations)
axConc = gca;
delta_axis  = 0.5;
axConc.YLim = [minConc-delta_axis maxConc+delta_axis];
axConc.XLim = [delta_axis length(plotStNames) + delta_axis];
set(axConc, 'XTick', 1:length(plotStNames), 'XTickLabel', plotStNames, 'XTickLabelRotation', LabelRotation, 'XAxisLocation', 'top')
axConc.Box = 'on';
axConc.YGrid = 'on';
axConc.XGrid = 'on';
axConc.YTick = minConc:1:maxConc;
axConc.GridLineStyle = '-';

% Draws rectangle for unfeasible concentration limits (above max and below min)
% e.g.: rectangle ('position',[xmin ymin xmin xmax])
concOutLim1 = rectangle('position',[0 -8 length(plotStNames)+1 2]);
set(concOutLim1, 'LineStyle', ':', 'LineWidth', 1, 'FaceColor', outerLimColor)

concOutLim2 = rectangle('position',[0 -2 length(plotStNames)+1 2]);
set(concOutLim2, 'LineStyle', ':', 'LineWidth', 1, 'FaceColor', outerLimColor)

%Draw Bottom X-axis and Right Y-Axis for pathway reactions
axConc_pos = axConc.Position;
ax2 = axes('Position', axConc_pos, 'XAxisLocation', 'bottom', 'YAxisLocation', 'right', 'Color', 'none');

%Sets characteristics of Bottom X-Axis and Right Y-Axis

ax2.YTick = '';
ax2.XLim = [0 length(reacNames)+1];
set(ax2, 'XTick', 1:length(reacNames), 'XTickLabel', '')

%% SUBPLOT 3: Proton Translocation and ATP plot
s3 = subplot(3,2,3); 
%Plots positive energy recovery sites
b1 = bar(protTranslocV' .* (posIndexes' == 1), 'FaceColor', posColor, 'BarWidth', barWidth);
hold on
%Plots reactions that need to be pumped in
b2 = bar(protTranslocV' .* (posIndexes' == 0), 'FaceColor', negColor, 'BarWidth', barWidth);
bATP = bar(n_ATPV * ratio_H_ATP, 'FaceColor', ATPColor, 'BarWidth', 0.5);

% hL_Energy = legend([b1(1),b2(1),bATP(1)],'H^+_{out}', 'H^+_{in}', 'ATP', 'Location', 'SouthEast','Orientation', 'Horizontal','box','off','FontSize',FontLegend);
%Set Transparency for bars
% b1.FaceAlpha = barTransparency;   b2.FaceAlpha = b1.FaceAlpha; bATP.FaceAlpha = b1.FaceAlpha; 
ylabel({'# of H^+';'Translocated'});
axProtRc = gca;
axProtRc.XLim  = ax2.XLim;
axProtRc.XTick = ax2.XTick;
axProtRc.XTickLabel = '';
axProtRc.YLim  = [min(nProt)-delta_axis max(nProt)+delta_axis];
axProtRc.YTick = min(nProt):deltaProt:max(nProt);
axProtRc.YGrid = 'on';
axProtRc.XGrid = 'on';
axProtRc.GridLineStyle = '-';

%% SUBPLOT 2: PLOT eC reoxidations 
%(Y Axis plotted in inverse fashion to make them coincidental with first plot)
s2 = subplot(3,2,2);

%Sets Transparency of bars
s1.FontSize = FontSize; ax2.FontSize = FontSize;
s2.FontSize = FontSize;
%Set Handle for eCarriers Left Y-Axis (Protons) and Top X-Axis
axProtons = gca;
delta_Xaxis = 0.5;
delta_Yaxis = 0;
axProtons.XLim = [delta_Xaxis length(plotRatios(1,:)) + delta_Xaxis];
axProtons.YLim = [nProt(1)-delta_Yaxis nProt(end)+delta_Yaxis];

set(axProtons, 'XTick', 1:length(plotRatios(1,:)), 'XTickLabel', eCarrierRatioLabel, 'YTickLabel', '', 'XAxisLocation', 'top')
axProtons.Box = 'off';      %sets off box to allow axis two to be used
axProtons.FontSize = FontSize;
% %If the carrier plotted is CoB, orientate the ratios
% if strcmp(eB_carrier, 'CoB-SH') 
%     axProtons.XTickLabelRotation = 30;
% else 
%     axProtons.XTickLabelRotation = 0;
% end
axProtons.XTickLabelRotation = LabelRotation;

%Draw next axis
axProtonspos = axProtons.Position;
axRatios = axes('Position', axProtonspos, 'XAxisLocation','bottom','YAxisLocation','right',...
    'Color','none');
hold on

%Plot Ratios
ratioPlot = plot(1:length(plotRatios(1,:)), log10(plotRatios)');
ylabel('log_{10}(Ratio_{e^- Carriers})')

%Set Marker Colours and Sizes
for i = 1:numRes
    ratioPlot(i).MarkerSize = 10;
    ratioPlot(i).Marker = markerList{i};
    ratioPlot(i).MarkerFaceColor = [0.8 0.8 0.8];
    ratioPlot(i).MarkerEdgeColor = [0   0   0];
    ratioPlot(i).LineStyle = 'none';
    ratioPlot(i).Color = [0 0 0];
end

axRatios.Box = 'off';
axRatios.YLim = [-4-1 4+1];  %Ratios limits
axRatios.XLim = axProtons.XLim;
axRatios.YGrid = 'on';
axRatios.XGrid = 'on';
axRatios.YTick = minRatio:deltaProt*2:maxRatio; 
axProtons.YTick = axRatios.YTick;
axProtons.YLim = axRatios.YLim;
set(axRatios, 'XTick', 1:length(plotRatios), 'XTickLabel', '')
% fix_xticklabels(axRatios);  %Write in two lines reaction names

% New rectangle outside of feasible Ratios eCarriers
concOutLim1 = rectangle('position',[0 -7 length(plotRatios(1,:))+1 3]);
set(concOutLim1, 'LineStyle', ':', 'LineWidth', 1, 'FaceColor', outerLimColor)

concOutLim2 = rectangle('position',[0 4 length(plotRatios(1,:))+1 3]);
set(concOutLim2, 'LineStyle', ':', 'LineWidth', 1, 'FaceColor', outerLimColor)

%% SUBPLOT 4: Subplot protons Translocated Ratios 
s4 = subplot(3,2,4);
%Plot Results
b3 = bar(bsxfun(@times, plotHTransloc, (plotHTransloc >= 0))', 'FaceColor', posColor, 'BarWidth', barWidth);
hold on 
b4 = bar(bsxfun(@times, plotHTransloc, (plotHTransloc < 0))', 'FaceColor', negColor, 'BarWidth', barWidth);

axProtRatios = gca; 
axProtRatios.YLim = axProtRc.YLim;
axProtRatios.YTick = axProtRc.YTick;
axProtRatios.YTickLabel = '';
axProtRatios.XTickLabel = '';
axProtRatios.YGrid = 'on';
axProtRatios.XGrid = 'on';
axProtRatios.GridLineStyle = '-';
axProtRatios.XTick = axProtons.XTick;
axProtRatios.XLim = axProtons.XLim;


%% SUBPLOT 5: Plot Dissipated Energy of each reaction step 
s5 = subplot(3,2,5);
b5 = bar(-DGrV', 'FaceColor', DGdissColor, 'BarWidth', barWidth);
axDG = gca;
axDG.XTick = ax2.XTick;
ylabel({' \DeltaG_{diss} ';' (kJ/mol)'})

if reacNameLabels == 1
    set(axDG, 'XTick', ax2.XTick, 'XLim', ax2.XLim, 'XTickLabel', reacNames)
else
    set(axDG, 'XTick', ax2.XTick, 'XLim', ax2.XLim, 'XTickLabel', reacPlotAbbrevNames)
end

maxDG = max(max(-DGrV));
DGrounder = roundn(maxDG,0);
axDG.YLim = [0 maxDG+1];
axDG.YTick = 0:round(DGrounder/2,1):DGrounder; 
axDG.YGrid = 'on';
axDG.FontSize = FontSize;
%Write names in two lines
if reacNameLabels == 1
    fix_xticklabels(axDG);
end


%% SUBPLOT 6: Plot Dissipated Energy of each reaction step 
s6 = subplot(3,2,6);
b6 = bar(-DGrV_eC', 'FaceColor', DGdissColor, 'BarWidth', barWidth);
axDG_eC = gca;
% ylabel({'\DeltaG_{diss} ';' (kJ/mol)'})
set(axDG_eC, 'XTick', axRatios.XTick, 'XLim', axRatios.XLim, 'XTickLabel', reacReoxAbbrevNames, ...
    'XTickLabelRotation', LabelRotation, 'YTickLabel', '')
axDG_eC.YLim = [0 max(max(-DGrV)) + 1];
axDG_eC.YTick = axDG.YTick;
axDG_eC.YGrid = 'on';
axDG_eC.XGrid = 'on';
axDG_eC.FontSize = FontSize;
%Write names in two lines
if reacNameLabels == 1
    fix_xticklabels(axDG_eC);
end

%Set Proper Positions for all figures and axis
s5.Position = [left_pos bottom_posDG, width_1, height_DG];
s2.Position = [left_pos + width_1 + plot_space, bottom_pos, width_2, height];
s1.Position = [left_pos, bottom_pos, width_1, height];
s6.Position = [left_pos + width_1 + plot_space, bottom_posDG,  width_2, height_DG];

s3.Position = [left_pos, bottom_protTransloc, width_1, height_HTransloc];
s4.Position = [left_pos + width_1 + plot_space, bottom_protTransloc, width_2, height_HTransloc];


axConc.Position =  [left_pos, bottom_pos, width_1, height];
ax2.Position =  [left_pos, bottom_pos, width_1, height];
axProtons.Position = [left_pos + width_1 + plot_space, bottom_pos, width_2, height];
axRatios.Position = [left_pos + width_1 + plot_space, bottom_pos, width_2, height];
s3.FontSize = FontSize;   axRatios.FontSize = FontSize;
axConc.TickLength = [0.005 0.015];
s3.TickLength = [0.005 0.015];


if reacLabelTicks == 0
%     ax2.TickLength  = [0 0];
    axDG.TickLength = [0 0];
%     axProtRc.TickLength = [0 0];
    ax2.XTick = '';
%     axDG.XTick = '';
    axProtRc.XTick = '';
else
    ax2.TickLength = axConc.TickLength;
    axDG.TickLength = axConc.TickLength;
    axProtRc.TickLength = axConc.TickLength;
end


%Add Extra Ticks
axProt_pos = axProtRc.Position;
ax3 = axes('Position',axProt_pos, 'XAxisLocation','bottom','YAxisLocation','right','Color','none');
ax3.XTick = axConc.XTick;
ax3.XLim  = axConc.XLim;
ax3.XTickLabel = '';
ax3.YTick = '';
ax3.TickLength =  [0.005 0.015];
ax3.XGrid = 'on';

axDG_pos = axDG.Position;
ax5 = axes('Position',axDG_pos, 'XAxisLocation','bottom','YAxisLocation','right','Color','none');
ax5.XTick = axConc.XTick;
ax5.XLim  = axConc.XLim;
ax5.XTickLabel = '';
ax5.YTick = '';
ax5.TickLength =  [0.005 0.015];
ax5.XGrid = 'on';

%Set Transparencies and Line Widths
bATP.FaceAlpha = barTransparency;   
bATP.LineWidth = barLineWidth;   

for i = 1:numRes
    
  b1(i).FaceAlpha = barTransparency;   
  b2(i).FaceAlpha = barTransparency;   
  b3(i).FaceAlpha = barTransparency;   
  b4(i).FaceAlpha = barTransparency; 
  b5(i).FaceAlpha = barTransparency; 
  
  b1(i).LineWidth = barLineWidth;   
  b2(i).LineWidth = barLineWidth;   
  b3(i).LineWidth = barLineWidth;   
  b4(i).LineWidth = barLineWidth; 
  b5(i).LineWidth = barLineWidth; 

end


