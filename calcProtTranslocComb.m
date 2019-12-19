function [protTransComb] = calcProtTranslocComb(protTranslocM, minProtTransloc, maxProtTransloc, n_H_ATP_SLP)

%Creates structure for every possible proton translocation at each step  
for i=1:length(protTranslocM(1,:))
Path.(char(strcat('Htrans_',num2str(i)))) = protTranslocM(1,i): protTranslocM(2,i);
end

%Collect all the names of the cells and add the string of the structure
pathNames = fieldnames(Path);
pathNames = strrep(pathNames,'H','Path.H');

%Creates a string and separates every cell name by commas
evalPaths = strjoin(pathNames, ',');

%Obtain all the possible combinations of proton translocatios
eval(strcat('protTransComb = allcomb(', evalPaths,');'))
% protTransComb = protPaths;

%Evaluate only the pathways that have a net proton translocation of bigger
%than -3 (Only one proton is produced via SLP, so a max of 3 can be
%translocated)
netProtTransloc = sum(protTransComb,2) + n_H_ATP_SLP;
%Proton Translocation combinations allowed from minimum to maximum protons
%allowed
protTransComb = protTransComb(and(netProtTransloc <= maxProtTransloc, netProtTransloc >= minProtTransloc),:);