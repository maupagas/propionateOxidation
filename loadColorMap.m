%This function creates a colormap for bars to be made sequentially

function [colorMap] = loadColorMap()

%********** EDIT COLORS FROM HERE *******************%

%Reds
redcolormap  = [130 0 0; 255 0 0; 255 100 100; 255 180 180];
redcolormap  = redcolormap/255;
%Blues
bluecolormap = [redcolormap(:,3), redcolormap(:,2), redcolormap(:,1)]; 
%Greens
greencolormap     = [30 102 54; 47 150 80; 100 210 135; 120 220 180]/255;
%Purples
purplecolormap    = [ 67  37  89; 136  75 181]/255;
%Turquoises
turquoisecolormap = [ 72 157 184; 138 193 210]/255;
%Oranges
orangecolormap    = [255 153   0; 255 189  91]/255; 

%********** EDIT COLORS UNTIL HERE *******************
%********* More colors can be added ******************


%Combine all colors previously defined
colorMap = [redcolormap; bluecolormap; greencolormap; purplecolormap; turquoisecolormap; orangecolormap];

