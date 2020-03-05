function [ConcV, prodConc, r1, r2, eBreax] = solve_eB(Param, StM, St, Reac, stoM, n_ATPV,  constConcM, protTranslocV, reacPath)

% id_St_eB = strcmp(St.StType, 'eB');
% St_eB = St.StNames(id_St_eB);
% 
% id_eB = find(strcmp(St_eB, Param.eBifurc));

id_Rc_eB = find(strcmp(Reac.(char(reacPath)).Labels_ord, 'eB'));

stoV = stoM(:, id_Rc_eB);
numProtTransloc = protTranslocV(id_Rc_eB);
n_ATP = n_ATPV(id_Rc_eB);

loopFlag = 0;    %No loop for this calculation

% %Find the reduced and the oxidised name of one of the two components of
% %the electron bifurcation
% eB_ox  = St_eB(id_eB);
% eB_red = St_eB(id_eB - 1);

%Find the reduced and the oxidised name of one of the ferredoxin
% id_Fdred = find(strcmp(St_eB, 'Fdred'));
% eC_ox  = St_eB(id_Fdred - 1);
% eC_red = St_eB(id_Fdred);
eBName      = strrep(Param.eB, '_', '-');
eB_ConsName = strrep(Param.eB_Cons, '_', '-');

%Find the position of all the components involved in the electron
%bifurcation
id_StM_eCred = (strcmp(St.StNames, Param.eB_eC));        %Ferredoxin reduced
id_StM_eCox  = (strcmp(St.StNames, Param.eB_eC_Cons));   %Ferredoxin oxidised
id_StM_eBred = (strcmp(St.StNames, eBName));           %CoM-SH
id_StM_eBox  = (strcmp(St.StNames, eB_ConsName));      %CoM-S-S-CoB

%Find reactions where bifurcations are involved
id_Rc_eC = (stoM(id_StM_eCred,:) ~= 0);
id_Rc_eB = (stoM(id_StM_eBred,:) ~= 0);
eBreax = id_Rc_eC + id_Rc_eB;

%Maximum ratio applicable 
r1 = Param.Max_Conc/Param.Min_Conc;

%Apply Maximum Ferredoxin Ratio
StM(id_StM_eCred) = StM(id_StM_eCox) * r1;

%Define constant concentrations matrix as to allow the reduced
%CoB/CoM-S-S-CoB ratio to be calculated
constConcV = ones(size(constConcM(:,1)));

%Make variable CoB-SH position
constConcV (id_StM_eBred) = 0;
constConcV = logical(constConcV);
varConcV = (constConcV == 0);

%Solve product concentration for CoB-SH
[ConcV, prodConc] = calcProdConc(Param, StM, stoV, n_ATP, constConcV, varConcV, numProtTransloc, loopFlag);            

StM(id_StM_eBred) = prodConc; 

%%%% Calculate ratio between reduced and oxidised e-bifurcator
r2 = StM(id_StM_eBred)/StM(id_StM_eBox);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check here how to cancel loop calculation if the ratio is below 1e-4   %%%
%%% (not enough concentration to run reaction forward)                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

