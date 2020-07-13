%This function creates a colormap for bars to be made sequentially

function [colorMap] = loadColorMap(plotAllFlag)

%********** EDIT COLORS FROM HERE *******************%
% create a default color map ranging from red to light pink
% redlength = 8;
% red_0 = [1, 0, 0];
% red_f = [255 180 180]/255;
% redcolormap = [linspace(red(1),pink(1),length)', linspace(red(2),pink(2),length)', linspace(red(3),pink(3),length)'];

%Reds
% redcolormap  = [130 0 0; 192.5 0 0; 255 0 0; 255 50 50; 255 100 100; 255 140 140; 255 180 180; 255 210 210];
% redcolormap  = redcolormap/255;
redlength = 8;
red_0 = [1, 0, 0];
red_f = [255 180 180]/255;
redcolormap = [linspace(red_0(1),red_f(1),redlength)', linspace(red_0(2),red_f(2),redlength)', linspace(red_0(3),red_f(3),redlength)'];

redcolormapL1 = [130 0 0]/255;
%Blues
% bluecolormap = [redcolormap(:,3), redcolormap(:,2), redcolormap(:,1)]; 

bluelength = 8;
blue_0 = [0, 0, 1];
blue_f = [180 180 255]/255;
bluecolormap = [linspace(blue_0(1),blue_f(1),bluelength)', linspace(blue_0(2),blue_f(2),bluelength)', linspace(blue_0(3),blue_f(3),bluelength)'];

bluecolormapL1 = [ 0 0 130]/255;
%Greens
% greencolormap     = [30 102 54; 47 150 80; 100 210 135; 120 220 180]/255;
greenlength = 8;
green_0 = [0.1, 0.41, 0.18];
green_f = [0.57, 1, 0.77];
greencolormap = [linspace(green_0(1), green_f(1),greenlength)', linspace(green_0(2),green_f(2),greenlength)', linspace(green_0(3),green_f(3),greenlength)'];

greencolormapL1 = [47 150 80]/255;
%Purples
% purplecolormap    = [ 67  37  89; 136  75 181]/255;
purplelength = 4;
% purple_0 = [67  37  89]/255;
% purple_f = [136  75 181]/255;
purple_0 = [0.39 0.11 0.45];
purple_f = [0.96 0.83 0.99];
purplecolormap = [linspace(purple_0(1), purple_f(1),purplelength)', linspace(purple_0(2),purple_f(2),purplelength)', linspace(purple_0(3),purple_f(3),purplelength)'];

purplecolormapL1 = [136  75 181]/255;

%Turquoises
turquoisecolormap = [ 72 157 184; 138 193 210]/255;

%Oranges
orangecolormap    = [255 153   0; 255 189  91]/255; 

%Grays
graycolormap    = [105,105,105; 128 128 128; 192 192 192; 220 220 220]/255;

%Browns
browncolomap    = [160 82 45; 210 105 30; 245 222 179; 255 248 220]/255;

%********** EDIT COLORS UNTIL HERE *******************
%********* More colors can be added ******************


%Combine all colors previously defined
if plotAllFlag == 1
    colorMap = [redcolormap; bluecolormap; greencolormap; purplecolormap; turquoisecolormap; orangecolormap; graycolormap; browncolomap];
else
    colorMap = [redcolormapL1; bluecolormapL1; greencolormapL1; purplecolormapL1; turquoisecolormap(2,:); orangecolormap(2,:); graycolormap(2,:); browncolomap(2,:)];    
end

