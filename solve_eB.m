function [ConcV, prodConc, r1, r2, eBreax] = solve_eB(Param, StM, St, Reac, stoM, n_ATPV,  constConcM, protTranslocV, reacPath)

id_St_eB = strcmp(St.StType, 'eB');
St_eB = St.StNames(id_St_eB);

id_eB = find(strcmp(St_eB, Param.eBifurc));

id_Rc_eB = find(strcmp(Reac.(char(reacPath)).Labels_ord, 'eB'));

stoV = stoM(:, id_Rc_eB);
numProtTransloc = protTranslocV(id_Rc_eB);
n_ATP = n_ATPV(id_Rc_eB);

loopFlag = 0;    %No loop for this calculation


%Find the reduced and the oxidised name of one of the two components of
%the electron bifurcation
eB_ox  = St_eB(id_eB);
eB_red = St_eB(id_eB - 1);

%Find the reduced and the oxidised name of one of the ferredoxin
id_Fdred = find(strcmp(St_eB, 'Fdred'));
eC_ox  = St_eB(id_Fdred - 1);
eC_red = St_eB(id_Fdred);

%Find the position of all the components involved in the electron
%bifurcation
id_StM_eCred = (strcmp(St.StNames, eC_red));
id_StM_eCox  = (strcmp(St.StNames, eC_ox));
id_StM_eBred = (strcmp(St.StNames, eB_red)); 
id_StM_eBox  = (strcmp(St.StNames, eB_ox));

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

%Solve product concentration for CoB-SH
[ConcV, prodConc, posConc] = calcProdConc(Param, StM, stoV, n_ATP, constConcV, numProtTransloc, loopFlag);            

StM(id_StM_eBred) = prodConc; 

%%%% Calculate ratio between reduced and oxidised e-bifurcator
r2 = StM(id_StM_eBred)/StM(id_StM_eBox);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check here how to cancel loop calculation if the ratio is below 1e-4 (not enough concentration to run reaction forward)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% switch loopFlag
%     
%     case 0 % There is no loop
%         %Solve sequentially the rest of the reaction
%          %Solve pathway with loop
%         for i = 1:length(protTranslocV)
%             stoV  = stoM(:,i);
%             n_ATP = n_ATPV(i);
%             constConcV = constConcM(:,i);
% 
%             %Calculate Ratio 1st Reaction
%             if strcmp(Reac.(char(reacPath)).Labels_ord(i),'Cr'),
%                 loopFlag = 2;
%             else
%                 loopFlag = 0;
%             end
%             %Calculate product concentration such as DGr = 0
%             [ConcV, prodConc, posConc] = calcProdConc(Param, StM, stoV, n_ATP, constConcV, protTranslocV(i), loopFlag);            
% 
%             %Store Results into structures
%             ConcResM(:,i)   = ConcV;
%             Prod_ConcV(i) = prodConc; 
%             Pos_ConcV(i)    = posConc;
% 
% 
%             %Updates the States Matrix with the calculated result
%             StM = ConcV;
%         end
%             
%     case 1   % Loop reactions, solve accordingly
%                 %Identify positions of the loop ****MOVE FROM HERE TO idProdsSubs()
% %                     [id, StM, stoLoopM, n_ATPLoopV, protTranslocLoopV, constConcLoopM] = identifyLoop(StM, Param, stoM, n_ATPV, protTranslocV, idCalcReax, St.id.H2);
% 
%                 %Collect how many loops are available in the pathway
%                 numLoops = idLoop.(char(reacPath)).numLoops;
% 
%                 calcLoops = 1;
%                 loopName = sprintf('Loop_%d', calcLoops);
% 
%                 counterRc = 1;
% 
%                 %Solve stoichiometry with loop
%                 while calcLoops <= numLoops || counterRc <= numReac
% 
%                     startLoop = idLoop.(char(loopName)).StartLoop;
%                     endLoop   = idLoop.(char(loopName)).EndLoop;
% 
%                     if counterRc < startLoop || counterRc > startLoop 
% 
%                         sprintf('Solve sequentially reaction %d', counterRc)
% 
%                         %Update matrices
%                         stoV  = stoM(:,counterRc);
%                         n_ATP = n_ATPV(counterRc);
%                         constConcV = constConcM(:,counterRc);
% 
%                         %Calculate product concentration such as DGr = 0
%                         [ConcV, prodConc, posConc] = calcProdConc(Param, StM, stoV, n_ATP, constConcV, protTranslocV(counterRc), counterRc, loopFlag);            
% 
%                         %Store Results into structures
%                         ConcResM(:,counterRc)   = ConcV;
%                         Prod_ConcV(j,counterRc) = prodConc; 
%                         Pos_ConcV(counterRc)    = posConc;
% 
%                         %Updates the States Matrix with the calculated result
%                         StM = ConcV;
%                         counterRc = counterRc + 1;
% 
%                     elseif counterRc == startLoop
%                             sprintf('Solve loop as a block reactions %d to %d', startLoop, endLoop)
% 
%                             %Pass up the information required for the loop
%                             LoopData = idLoop.(char(reacPath)).(char(loopName));
%                             %Solve Loop
%                             [ConcV, LoopResults, solutionFlag] = solveLoop(StM, Reac, Param, LoopData);
% 
%                             %Store Results into structures
%                             ConcResM(:, counterRc : id.EndLoop)   = LoopResults.ConcV;
%                             Prod_ConcV(j, counterRc : id.EndLoop) = LoopResults.Prod_Conc; 
%                             Pos_ConcV(:, counterRc : id.EndLoop)  = LoopResults.Pos_Conc;
% 
%                             %Updates the States Matrix with the calculated result
%                             StM = ConcV;                              
% 
%                             %Update Loop counter                                
%                             counterRc = endLoop + 1;
%                             % Update number of loops calculated counter
%                             calcLoops = calcLoops + 1;
%                             if calcLoops <= numLoops,
%                                 loopName = sprintf('Loop_%d', calcLoops);    
%                             end
%                     end
%                 end
% 
%                             sprintf('Whole stoichiometry with loop has been solved')
%         
%        
%         
%         
%         
%         
%         
% end
%     
% 
%     [DGrV] = calcEnergetics(StM, Param, Reac, stoM, n_ATPV, protTranslocV);
% 
%     if all(DGrV > 1e-4) == 1,
%         solutionFlag = 0;
%     elseif all(DGrV < 1e-4),
%         solutionFlag = 1;
%     end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

%calcProdConc

