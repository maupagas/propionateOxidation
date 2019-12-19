function [] = plotTable(St, numResults, FinalResultsSorted, netATPSorted)

%% WRITE TABLE WITH RESULTS

    %Labels and Inputs results to plot Table ,
    cnames  = St.combinationListNames;
    colRatio = find(strcmp(cnames,'ratio_H_ATP'));
    
    % Write cnames as HTML code for table usage
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
    
    rnamesTable = cell(length(numResults),1);
    % rnamesLegend = rnamesTable;
    
    % rnamesLegend1  = strcat('P_{H2} = ', strData(:,col_H2),'ppm');
    % rnamesLegend2 = strcat('Net ATP = ', strData(:,end));
    
    % rnamesLegend = sprintf('%s\n%s', char(rnamesLegend1), char(rnamesLegend2));
    
    % ResData =  PrintResults.FeasComb.(char(reacPath)).Inputs(numResults,:);
    % netATP  =  PrintResults.FeasComb.(char(reacPath)).netATPV(numResults,:);
    
    ResData = [FinalResultsSorted(numResults,:), netATPSorted(numResults,:)];
    
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
                P_H2 = ResData(i,j) /(exp(-DG_h2/(Rth * (T +35)))) * 1e5;    %Pa
                strData{i,j} = sprintf('%.2f', P_H2);
            end
        end
    end
    
    
    for i = 1:length(numResults)
        rnamesTable{i}  = sprintf('<html><tr align=center><col width=80><font size=+1>Config<sub>%d</sub>',i);
        rnamesLegend1 = sprintf('P_{H2} = %s Pa', char(strData(i,col_H2)));
        rnamesLegend2 = sprintf('Net ATP = %s', char(strData(i,end)));
        
        rnamesLegend{i} = sprintf('%s\n%s', char(rnamesLegend1), char(rnamesLegend2));
        
    end

%% ************************************************
%************ TABLE DIMENSIONS ******************
%************************************************
t_left   = 0.07;
t_bottom = 0.82;
t_width  = 0.095 * length(ResData(1,:));

if length(numResults) >= 3
    t_height = 0.06 * length(numResults);
elseif length(numResults) == 2
    t_height = 0.07 * length(numResults);
elseif length(numResults) == 1
    t_height = 0.1 * length(numResults);
end

HeaderFontSize = 1;   HeaderWidth = 10;   columnWidth = 115;
strCenteredData = strcat(sprintf('<html><tr align=center><td width=%d>', columnWidth), strData(:,:));
cnames = strcat(sprintf('<html><font size=+%d>',HeaderFontSize), cnames(:));
% rnamesHeader = strcat(sprintf('<html><font size=+%d>',HeaderFontSize), rnamesTable(:));
rnamesHeader = rnamesTable;
r_H_ATP_pos   = strcmp(cnames, 'ratio_H_ATP');
ratio_H_ATP   = ResData(1,r_H_ATP_pos);

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
