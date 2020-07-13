function []=printFile(outputFileName, printPath, paperWidth, paperHeight)

%This function is built to save the plots as MATLAB Figure .fig, as PDF and as image .png
 
%Set the units of the paper
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);
set(gcf, 'PaperSize', [paperWidth paperHeight]); % dimension on x axis and y axis resp.

%Save Figure, image as pdf and png for publication purposes
savefig(gcf, sprintf('%s%s.fig',       printPath, outputFileName))                  %.fig file
print(gcf,'-dpdf', sprintf('%s%s.pdf', printPath, outputFileName))                  %.pdf file  
print(gcf,sprintf('%s%s.png',          printPath, outputFileName),'-dpng','-r600')  %.png file  

%Alternative
% exportgraphics(gcf, sprintf('%s%s.fig', printPath, outputFileName), 'Resolution',300)

