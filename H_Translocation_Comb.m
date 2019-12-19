function [ProtTransComb] = calcProtTranslocComb(ProtonTranslocM, minProtTransloc, maxProtTransloc)

%Creates structure for every possible proton translocation at each step  
for i=1:length(ProtonTranslocM(1,:))
Path.(char(strcat('Htrans_',num2str(i)))) = ProtonTranslocM(1,i):ProtonTranslocM(2,i);
end

%Collect all the names of the cells and add the string of the structure
PathNames = fieldnames(Path);
PathNames = strrep(PathNames,'H','Path.H');

%Creates a string and separates every cell name by commas
EvalPaths = strjoin(PathNames, ',');

%Obtain all the possible combinations of proton translocatios
eval(strcat('H_Paths = allcomb(', EvalPaths,');'))
ProtTransComb = H_Paths;

%Evaluate only the pathways that have a net proton translocation of bigger
%than -3 (Only one proton is produced via SLP, so a max of 3 can be
%translocated)
NetHTransloc = sum(H_Paths,2);
ProtTransComb = ProtTransComb(NetHTransloc > -3,:);