%Define colormap to use
redcolormap  = [130 0 0; 255 0 0; 255 100 100; 255 180 180];
% bluecolormap = [redcolormap(:,3), redcolormap(:,2), redcolormap(:,1)]; 
redcolormap  = redcolormap/255;
bluecolormap = [redcolormap(:,3), redcolormap(:,2), redcolormap(:,1)]; 
% greencolormap = [redcolormap(:,3), redcolormap(:,1), redcolormap(:,2)]; 
greencolormap     = [30 102 54; 47 150 80; 100 210 135; 120 220 180]/255;
purplecolormap    = [ 67  37  89; 136  75 181]/255;
turquoisecolormap = [ 72 157 184; 138 193 210]/255;
orangecolormap    = [255 153   0; 255 189  91]/255; 

pathwaymap = [redcolormap; bluecolormap; greencolormap; purplecolormap; turquoisecolormap; orangecolormap];

% colormap(pathwaymap);
reacPathway = Reac.Pathway(Param.firstPath2Eval:Param.lastPath2Eval);

%Define starting point of the barplot
minProtTransloc = -1;
sensTitle = 'Sensitivity analysis';
%Variables of the name to plot
outputVarNames = {'\DeltaG_{ATP}  (kJ/mol)', 'Num H^+ / ATP', 'CoA-SH  (M)', 'pH_{in}', 'P_{H2} (ppm)', 'Dissolved [CO_2]  (M)', 'Temperature (^oC)', 'pH_{out}'};

%Convert units
T = 273.15;
Rth = 8.314e-3;
DG_h2 = 18.3;
DG_ATP = abs(Param.refComb(1));


%Label for Net Proton Translocation Axis
deltaPlot = 0;
minYval = 0;
maxYval =  3 + 1/3;
minYPlotVal = minYval/3 - deltaPlot;
maxYPlotVal = (DG_ATP+5)/DG_ATP + deltaPlot;

yticknumvals = (minYval:1/3:maxYval);
denomTicks = 3;
denomRatio = 3;

%Provide rational values for Y axis
[ytickvals] = rationalTickLabels(yticknumvals, denomTicks, 1);



% %Plot Results
for resNum = 1:2
    
%     if resNum == 1
%         load ProOxidH2toPlot.mat
%     else
%         clearvars Param PrintResults St Reac
%         load H2mtoPlot.mat
%     end
    var2Extract = 'H2';
    outputVar = var2Extract;
    
   if length(Output.varValues.(char(outputVar))) > 1
        plotResults = Output.plotResults.(char(outputVar));
        
        xTickLabel = Output.varValues.(char(outputVar));
        
        if strcmp(outputVar,'H2') == 1
            P_H2 = Output.varValues.H2 /(exp(-DG_h2/(Rth * (T +35)))) * 1e6;
            for n = 1:length(P_H2)
                P_H2_Label{n} = sprintf('%.1f',P_H2(n));
            end
            xTickLabel = P_H2_Label;
        end
        
        %If there is more than 4 plots, open a new figure
        if resNum == 1 
            figure('units','normalized','outerposition',[0 0 1 1]);
        end
        
        subplotNum = resNum;
        
        %Establishes 2 plots per figure
        s = subplot(2,1,subplotNum);
        b = bar(plotResults,'BaseValue', minProtTransloc);
        ax1 = gca;
        ax1.YGrid = 'on';
        
        % For MATLAB 2018, Apply each color to each bar
        for y = Param.firstPath2Eval:Param.lastPath2Eval 
            y =  y - Param.firstPath2Eval + 1;
            b(y).FaceColor = pathwaymap(y,:);
            b(y).EdgeColor = 'k';
        end
        
        %Ticks the labels for each variable in the y-axis
%         if strcmp(outputVar, 'ratio_H_ATP') == 1
%             [ratio_H_Label] = rationalTickLabels(Output.varValues.ratio_H_ATP, denomRatio, 0);
%             set(gca, 'XTick', 1:length(plotResults(:,1)), 'XTickLabel', ratio_H_Label);
%         else
%
%         end
        %     set(gca, 'YTickLabel', Output.varValues.(char(outputVar)));
        %     set(gca, 'Ydir','reverse')
        set(ax1, 'XTick', 1:length(plotResults(:,1)), 'XTickLabel', xTickLabel);
        set(ax1, 'YTick', yticknumvals, 'YTickLabel', ytickvals);
        set(ax1, 'Ylim',[minYPlotVal maxYPlotVal]);
        
        %Define measures of plots
        plotWidth  = 0.80;
        plotHeight = 0.37;
        xPos = 0.10;
        yPos = 0.60;
        gapHeight = 0;
        gapWidth  = 0.005;
        set(gca,'FontSize',14)
        
        ax1.Position = [0.15 0.65 0.75 0.3];
%         ax2.Position = [0.15 0.65 0.75 0.3];
        
        %Define position of subplots
        if subplotNum == 1
            Position = [xPos, yPos, plotWidth, plotHeight];
            %         title(sensTitle, 'Position', [2.05 0.05 0]);
            ylabel({'Net ATP_{harvested}','(Y_{ATP/pro})'})
            hL = legend(reacPathway,'NumColumns',9);
            hL.Position = [0.3 0.05 0.45 0.05];

        elseif subplotNum == 2
            Position = [xPos, yPos-plotHeight-gapHeight, plotWidth, plotHeight];
            ylabel({'Net ATP_{harvested}','(Y_{ATP/H2})'})
        end
        set(ax1, 'Position',Position);

        %Collect values of DG
        DGr_vals = zeros(1,length(xTickLabel));
        for l = 1:length(xTickLabel)
            varDGr_Name = strcat(outputVar,'_',num2str(l));
            DGr_vals(l) = Output.DGr.(char(varDGr_Name));
        end
        
        %Prepare to plot the DGr on the secondary axis
        ax1.GridColor = 'k';
        ax1.GridLineStyle = ':';
        ax1.GridAlpha = 0.7;
        set(ax1,'box','off');
        if resNum == 2
            ax1.XLim = plot1_ax1XLim;
        end
          
        % close al
        
%         if resNum ~= 1
ax2 = axes('Position',ax1.Position, 'XAxisLocation','top','YAxisLocation','right','Color','none');

ax2.YLim = ax1.YLim * DG_ATP;
ax2.XLim = ax1.XLim;
set(ax2,'XTick' , 1:length(xTickLabel))
set(ax2,'XTickLabel' , '')
set(ax2,'YTick', ax1.YTick * DG_ATP)
% ax1.XLimMode = 'manual';
set(ax2,'YTickLabel', round(yticknumvals * DG_ATP,1))
ax2.YLabel.String = {'-\DeltaG_{Cat, available}' ;'(kJ/mol)'};
plot1_ax1XLim = ax1.XLim;

ax1.XLabel.String='P_{H2} (ppm)';
hold on
%     hDG = plot(-DGr_vals,'ok','MarkerFaceColor',[.85 .85 .85], 'LineWidth', 1, 'MarkerSize', 12);
sizeErrBar = 0.4;
errorbarSize = ones(1, length(xTickLabel)) * sizeErrBar;
hDG_error = herrorbar(1:length(xTickLabel), -DGr_vals, errorbarSize, errorbarSize, 'ok');

set(hDG_error, 'MarkerSize',  2, 'MarkerFaceColor','k', 'LineWidth', 1);
set(ax2,'FontSize',14)
if resNum == 1, legend('\DeltaG','Location', 'NorthEast','Color','w'); end
            %     hDG_error.MarkerFaceColor = [.85 .85 .85];
            %     hDG.LineWidth = 1.5;
            
            
%         else
%              ax2 = axes('Position',ax1.Position, 'XAxisLocation','top','YAxisLocation','right','Color','none');
%              ax2.YTick = '';
%              ax2.XTick = 1:length(xTickLabel);
%              ax2.XLim = ax1.XLim;
%              ax2.XTickLabel = '';            
%         end
        
        %     set(hleg,'Position',[0.28 -0.4841 0.450 0.5883])
        
   
    end
end


    
