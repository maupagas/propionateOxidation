%This script is written to plot the pathways of level two (different electron carriers selected)
%with highest ATP Yields. 

%Print Results Properties
close all
printToFile = 1;
printPath = strcat(pwd, '\Results\');
% Print all figures
paperWidth   = 24;
paperHeight = 16;
%Option to plot simple bars for Level 1 PAthway (flag = 0) or plot all
%(flag = 1)
plotAllFlag = 0;
reacPathL1 = Reac.reacPathL1;
% outputFileName = 'Scenario';

% Load Reac Pathway names
reacPath = Reac.Pathway(Param.firstPath2Eval:Param.lastPath2Eval);

%Define starting point of the barplot
minProtTransloc = -1;
sensTitle = 'Sensitivity analysis';
% Variables of the name to plot
outputVarNames = {'\DeltaG_{ATP}  (kJ/mol)', 'Ratio H^+ / ATP', '[CoA-SH]  (mM)', 'pH_{in}', 'P_{H2} (Pa)', '[CO_2]_{Dissolved} (mM)', 'Temperature (^oC)', 'pH_{out}'};

%Write the output yield units (mol ATP/mol eD)
eDvar = St.StNames(St.eD == 1);
eDvar = strrep(eDvar, '_out', '^-');
unitsYAxis = sprintf('(%s %s)', 'mol ATP/mol', char(eDvar));

% Values used to convert units
T_Kelvin = 273.15;   %For temperature to Celsius
Rth = 8.314e-3;      %For conversions of H2 to Pa
id_T = strcmp(Param.combNames, 'T');
id_ATP = strcmp(Param.combNames, 'DG_ATP');
DG_h2 = (Param.H0fM(St.id.H2) + (Param.refComb(id_T) / 298.15)*(Param.G0fM(St.id.H2) - Param.H0fM(St.id.H2)));   %Watch out when this is calculated for ALL POSSIBLE COMBINATIONS
DG_ATP = abs(Param.refComb(id_ATP));
fontSize = 12;
numPaths = length(Reac.Path(:,1));

%% Label for Net Proton Translocation Axis
deltaPlot = 0;
minYval = 0;
minYPlotVal = minYval/3 - deltaPlot;

maxYval =  Output.maxATP + 2/3;
maxYPlotVal = -Output.maxDGr/DG_ATP + deltaPlot;

yticknumvals = (minYval:1/3:maxYval);
denomTicks = 3;
denomRatio = 3;



%Define size of error bars as a function of the number of plots 
if plotAllFlag == 1
    sizeErrBar = 0.02;
else
    sizeErrBar = 0.06;
end

% Load colors for the plots
[pathwayMap] = loadColorMap(plotAllFlag);

% Provide rational values for Y axis
[ytickvals] = rationalTickLabels(yticknumvals, denomTicks, 1);

%% Plot Results for each variable
for i = 1:length(Param.combNames)
    
    %Select variable to plot 
    outputVar = Param.combNames(i);
    
    %If the variable is 
    if length(Output.varValues.(char(outputVar))) > 1
        
        if plotAllFlag == 1
            plotResults = Output.plotResults.(char(outputVar));
        else 
            plotResults = Output.plotResultsL1.(char(outputVar));
        end
            
        xTickLabel = Output.varValues.(char(outputVar));
        
        %Convert labels to desired axis values (e.g. H2 in ppm, T in C, etc)
        if strcmp(outputVar,'H2') == 1
            P_H2 = Output.varValues.H2 /(exp(-DG_h2/(Rth * (T_Kelvin +35)))) * 1e5;  %Pa
            P_H2_Label = cell(1, length(P_H2));
            for n = 1:length(P_H2)
                if P_H2(n)
                    P_H2_Label{n} = sprintf('%.2f',P_H2(n));
                else
                    P_H2_Label{n} = sprintf('%.1f',P_H2(n));
                end
            end
            xTickLabel = P_H2_Label;
        elseif strcmp(outputVar,'T') == 1
            T_Label = Output.varValues.T - T_Kelvin;
            xTickLabel = T_Label;
        elseif strcmp(outputVar,'CO2') == 1 || strcmp(outputVar,'CoA_SH') == 1
            xLabel = cell(1,length(xTickLabel));
            for n = 1:length(Output.varValues.(char(outputVar)))
                xLabel{n} = sprintf('%.1d',Output.varValues.(char(outputVar))(n)*1000);
            end
            xLabel = strrep(xLabel, '.','');
            xLabel = strrep(xLabel, 'e-0','^-^');
            xTickLabel = xLabel;
        end
        %***************** Labelling ends here ****************************
        
                
        %If there is more than 4 plots, open a new figure
        if i == 1 || i == length(Param.combNames)/2 + 1
            figure('units','normalized','outerposition',[0 0 1 1]);
        end
        
        subplotNum = mod(i, 4);
        if subplotNum == 0
            subplotNum = 4;
        end
        
        %Establishes 2 plots per figure
        subplot(2,2,subplotNum)
        
        b = bar(plotResults,'BaseValue', minProtTransloc);
        ax1 = gca;
        ax1.YGrid = 'on';
        
        % From MATLAB 2018b ONWARDS, Apply each color to each bar
        for y = 1:length(plotResults(1,:))
            b(y).FaceColor = pathwayMap(y,:);
            b(y).LineWidth = 0.1;
            b(y).EdgeColor = [0 0 0];
%             if strcmp(outputVar, 'Ratio_H_ATP') == 1
%                 b(y).EdgeColor = [0 0 0];
%             end            
        end
        
        %Ticks the labels for each variable in the y-axis
        if strcmp(outputVar, 'Ratio_H_ATP') == 1
            [ratio_H_Label] = rationalTickLabels(Output.varValues.Ratio_H_ATP, denomRatio, 0);
            set(gca, 'XTick', 1:length(plotResults(:,1)), 'XTickLabel', ratio_H_Label);
        else
            set(gca, 'XTick', 1:length(plotResults(:,1)), 'XTickLabel', xTickLabel);
        end
        %     set(gca, 'YTickLabel', Output.varValues.(char(outputVar)));
        %     set(gca, 'Ydir','reverse')
        set(gca, 'YTick', yticknumvals, 'YTickLabel', ytickvals);
        set(gca,'Ylim',[minYPlotVal maxYval]);
        
        %Define measures of plots
        plotWidth  = 0.35;
        plotHeight = 0.30;
        xPos = 0.15;
        yPos = 0.65;
        nCols = ceil(length(Reac.Path(:,1))/2);
        gapHeight = 0.11;
        gapWidth  = 0.005;
        set(gca,'FontSize',fontSize)
        
        %Define position of subplots
        if subplotNum == 1
            Position = [xPos, yPos, plotWidth, plotHeight];
            %         title(sensTitle, 'Position', [2.05 0.05 0]);
            ylabel({'Y_{ATP}',unitsYAxis})
            if plotAllFlag == 1
                nCols = ceil(length(Reac.Path(:,1))/4);
                hL = legend(reacPath,'NumColumns',nCols, 'EdgeColor', 'none');
            else
                hL = legend(reacPathL1,'NumColumns',nCols, 'EdgeColor', 'none');                
            end
            
            hL.Position = [0.3 0.05 0.45 0.05];
        elseif subplotNum == 2
            Position = [xPos+plotWidth+gapWidth, yPos, plotWidth, plotHeight];
            set(gca,'YTickLabel',[]);
        elseif subplotNum == 3
            Position = [xPos, yPos-plotHeight-gapHeight, plotWidth, plotHeight];
            ylabel({'Y_{ATP}',unitsYAxis})
        else
            Position = [xPos+plotWidth+gapWidth, yPos-plotHeight-gapHeight, plotWidth, plotHeight];
            set(gca,'YTickLabel',[]);
            %         legend(Reac.Pathway,'Location', 'BestOutside');
        end
        set(gca, 'Position',Position);
        textLabel = strcat(char(96+i),')');
        hText = text(0.03, 0.93 , textLabel, 'Units', 'normalized', 'FontSize', fontSize); 
%         hText.Position = [0.03 0.93 0];   % Defined position of the legend by hand

        %**************** POSITIONS DEFINITION ENDS HERE *******************

        %Collect values of DG
        if plotAllFlag == 1
            DGr_vals = zeros(length(xTickLabel),numPaths);
            for l = 1:length(xTickLabel)
                varDGr_Name = strcat(outputVar,'_',num2str(l));
                DGr_vals(l,:) = Output.DGr.(char(varDGr_Name));
            end
        else
            DGr_vals = Output.DGrResultsL1.(char(outputVar));
        end
            

        
        %Plot Gibbs Free Energies on the secondary axis
        if i ~= 1
            ax2 = axes('Position',ax1.Position, 'XAxisLocation','top','YAxisLocation','right','Color','none');
            
            ax2.YLim = ax1.YLim * DG_ATP;
            ax2.XLim = ax1.XLim;
            set(ax2,'XTick' , 1:length(xTickLabel))
            set(ax2,'XTickLabel' , '')
            set(ax2,'YTick', ax1.YTick * DG_ATP)
            ax1.XLimMode = 'manual';
            
            if mod(i,2) ~= 0
                set(ax2,'YTickLabel' , '')
                ax2.YLabel.String = '';
            else
                set(ax2,'YTickLabel', round(yticknumvals * DG_ATP,1))  %Put the tick values on the secondary y-axis
                ax2.YLabel.String = {'-\DeltaG_{Cat} (kJ/mol)'};
            end            
            hold on
            
%             sizeErrBar = 1/20; %Error bars as the horizontal plot lines
%             errorbarSize = ones(1, length(xTickLabel)) * sizeErrBar;
%             hDG_error = herrorbar(0.5:1:length(xTickLabel)-0.5, -DGr_vals, errorbarSize, errorbarSize, 'ok');
%           
%             set(hDG_error, 'MarkerSize',  0.2, 'MarkerFaceColor','k', 'LineWidth', 1);
            errorVals = -DGr_vals';
            %%For MATLAB 2019b or later releases
            % Calculate the number of bars in each group
            nbars = size(plotResults, 2);
            % Get the x coordinate of the bars
            xbars = zeros(nbars, length(xTickLabel));
            for n = 1:nbars
%                 x = [x ; b(n).XEndPoints];
                xbars(n,:) = b(n).XEndPoints;
            end


%             set(hDG_error2, 'MarkerSize',  0.2, 'MarkerFaceColor','r', 'LineWidth', 1, 'LineStyle', '-.');


            set(ax2,'FontSize',fontSize)
            
        else
            ax2 = axes('Position',ax1.Position, 'XAxisLocation','top','YAxisLocation','right','Color','none');
            ax2.YTick = '';
            ax2.XTick = 1:length(xTickLabel);
            ax2.XLim = ax1.XLim;
            ax2.XTickLabel = '';
            
%             nbars = size(plotResults, 2);
%             % Get the x coordinate of the bars
%             x = [];
%             for n = 1:nbars
%                 x = [x ; b(n).XEndPoints];
%             end
%             % Plot the errorbars
%             errorbar(x',errorVals',0.001*errorVals','k','linestyle','none')'
            
        end
        
        %Put Tick in between labels
        ax1Lim = ax1.XLim;
        deltaVars = (ax1Lim(2) - ax1Lim(1) ) / length(xTickLabel) /2;
        ax2.XLim = ax1.XLim - deltaVars;
        ax2.XGrid = 'on';
        ax2.XAxisLocation = 'bottom';
        ax1.Box = 'on';
        ax1.TickLength = [0 0];      
        ax1.XLabel.String=(outputVarNames(i));
        
        % Plot the errorbars
        if i ~= 1
%             hErr = errorbar(xbars'-deltaVars,errorVals',0.001*errorVals','k','linestyle','none');
%             set(hErr , 'MarkerSize',  0.2, 'MarkerFaceColor','k', 'LineWidth', 1);
            
%                 errorbar(xbars-deltaVars,errorVals,0*errorVals+0.01,0*errorVals+0.01);
%                 set(hErr2, 'MarkerSize',  0.2, 'MarkerFaceColor','k', 'LineWidth', 1);
            for m = 1:numel(plotResults)
                hErr2 = herrorbar(xbars(m)-deltaVars,errorVals(m),0*errorVals(m)+sizeErrBar,0*errorVals(m)+sizeErrBar, 'ok');
                set(hErr2, 'MarkerSize',  0.2, 'MarkerFaceColor','k', 'LineWidth', 1);
            end

            %Set legend of DG 
            if mod(i,4) == 0
                legend('\DeltaG','Color','w'); 
            end
        end
    end
end

if printToFile == 1
    
    if plotAllFlag == 1
        outputFileName = 'Y_ATP_';
    else
        outputFileName = 'Y_ATP_L1_';        
    end
    % Prepare print of the figures as PDF and JPG
    for j = 1:2
        figure(j)
        outputName = sprintf('%s%s', outputFileName, num2str(j));
        printFile(outputName, printPath, paperWidth, paperHeight)
    end
end

clearvars -except Reac Param idLoop St Results PrintResults Combination Output
