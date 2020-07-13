% A branched reaction is defined by the following reaction:
% $C_{in}  <----> C_{out} + C_1$
%
% $C_{in}  <----> C_{out} + C_1$
% $\Delta G_R = 0, C_{out} = 10^{-6} M$ for calculating C_{1,max} and minimum ratio

function [StM, LoopResults, solutionFlag] = solveBranched(StM, Param, LoopData, protTranslocLoopV, nProtComb)

%Preload States
DG0ft    = Param.DG0ft;        Rth      = Param.Rth;          T  = Param.T;        
DG_Prot  = Param.DG_Prot;      DG_ATP   = Param.DG_ATP;
Max_Conc = Param.Max_Conc;     Min_Conc = Param.Min_Conc;

%Load matrices
stoLoopM       = LoopData.stoLoopM;      
n_ATPLoopV     = LoopData.n_ATPLoopV;    
constConcLoopM = LoopData.constConcLoopM;
varConcLoopM   = LoopData.varConcLoopM;
posConcLoopV   = LoopData.posConcLoopV;
% protTranslocLoopV = LoopData.protTranslocLoopV(nProtComb,:);

%Identifiers of elements of the loop
id_C1 = LoopData.id_C1;              id_Cout = LoopData.id_Cout;
% id_Cn = LoopData.id_Cn;
% startLoop = LoopData.StartLoop;
%
nonC1V   = [1:id_C1-1   id_C1+1:length(DG0ft)];
nonCoutV = [1:id_Cout-1 id_Cout+1:length(DG0ft)];
flagCalcR1 = 0;   % This is 0 if R1 has not been previously calculated
sum_DG_loop = 0;   %Initialiser for this value

%INITIALIZE FIRST WHILE LOOP
concCheck  = 1;         DGCheck = 1; 
Cout_Check = 0;
C1 = Max_Conc;        
% Cout = Min_Conc;
DGTol = 1e-6;
loopFlag = 0; % No high concentrations are allowed when solving branched reactions

%Concentration of Cout is the lowest possible to give 'all the energy
%available' to the Loop
StM(id_C1) = C1;

%Preallocate LoopResults
LoopRes_ConcV     = zeros(length(StM), length(stoLoopM(1,:)));
LoopRes_Prod_Conc = zeros(1,length(stoLoopM(1,:)));
LoopRes_Pos_Conc  = LoopRes_Prod_Conc;

%% Calculate concentration of first component
while concCheck == 1 || DGCheck == 1 || Cout_Check == 0
    %         if ratioNew == 0,  sprintf('First calculation is being performed. \n'), else sprintf('Ratio is smaller than expected ratio. \n', sum(DGrVLoop)), end
    c_i_Loop = 0 * stoLoopM(1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%  CALCULATE LOOP CONCENTRATIONS  %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:length(stoLoopM(1,:))
        stoV = stoLoopM(:,i);
        ObjDG = 0 - protTranslocLoopV(i) * DG_Prot;
        n_ATP_Loop = n_ATPLoopV(i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % A. SOLVE FIRST REACTION BY PROVIDING A VALUE FOR C1 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i == 1
            
            %If it is the first calculation, no new ratio is provided
            if flagCalcR1 == 0
                              
                %Calculate Cout assuming a C1 as maximum concentration
                %(most favorable scenario for long pathway)
                u_i = DG0ft(nonCoutV) + Rth * T * log(StM(nonCoutV));
                sum_u_i = sum(u_i .* stoV(nonCoutV)) + ObjDG - DG_ATP * n_ATP_Loop;
                
                %Calculates the Product concentration (and the initial ratio) for a DG = 0
                Cout = exp((-sum_u_i/stoV(id_Cout) - DG0ft(id_Cout)) / (Rth * T));
                
                %If C1 is bigger than the maximum allowed concentration,
                %the reaction loop has to dissipate energy (If DG = 0 with
                %C1 > Max_Conc, it will be also favourable if C1 <
                %Max_Conc, as C1 is a product of the loop reaction
                if Cout > Max_Conc 
                    Cout = Max_Conc; 
                elseif Cout < Min_Conc
                    break
                end
                
                % Reaction 1 has been calculated, update the flag
                flagCalcR1 = 1;
            else    %For the rest of calculations in the loop, C1 concentration is provided
                Cout = Cout_new;
                
                %Calculate C1 by using the new value of Cout
                u_i = DG0ft(nonC1V) + Rth * T * log(StM(nonC1V));
                sum_u_i = sum(u_i .* stoV(nonC1V)) + ObjDG - DG_ATP * n_ATP_Loop;
                
                %Calculates the Product concentration (and the initial ratio) for a DG = 0
                C1_new = exp((-sum_u_i/stoV(id_C1) - DG0ft(id_C1)) / (Rth * T));
                
                %Update states vector
                StM(id_C1) = C1_new;
            end
            
            %Update states vector
            StM(id_Cout) = Cout;
            prodConc = Cout;
            LoopRes_Prod_Conc(i) = prodConc;
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % B. CALCULATE SEQUENTIALLY REST OF THE PRODUCT CONCENTRATIONS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            
%             %To allow the concentration of Cn to go over the physiological
%             %limits only when the ratio has not been determined that it
%             %will be bigger
%             if i == length(stoLoopM(1,:)) && ratioFlag < 2
%                 loopFlag = 1;
%             else
%                 loopFlag = 0;
%             end
%             
            %Calculate product concentration of the rest of the components
            
            constConcV = constConcLoopM(:,i);
            varConcV = varConcLoopM(:,i);
            posConc = posConcLoopV(i);

            [ConcV, prodConc] = calcProdConc(Param, StM, stoV, n_ATP_Loop, constConcV,  varConcV, protTranslocLoopV(i), loopFlag);
                        
            %Store Results into structures
            LoopRes_ConcV(:,i)   = ConcV;
            LoopRes_Prod_Conc(i) = prodConc;
            LoopRes_Pos_Conc(i)  = posConc;
                      
            %Evaluate the branched reaction once we reach the last reaction
            if i == length(stoLoopM(1,:))
                Cout_new = prodConc;
                
                %If the new Cout is lower than the initial value, the
                %calculation ends. If it is bigger, we need to re-calculate
                if Cout_new < Min_Conc
                    continue
                elseif Cout_new <= Cout && Cout_new > Min_Conc
                    StM = ConcV;
                    solutionFlag = 1;
                    Cout_Check = 1;
                else
                    Cout = Cout_new;
                    solutionFlag = 0;
                end
            else
                StM = ConcV;
            end
            
            %Prints the calculation for one component in the loop
            %sprintf('Concentration of %s is %.2d M.\n', char(StNames(posConc)), prodConc)
            
        end
        %Collect product concentrations in a vector
        if i == 1
            c_i_Loop(i) = C1;
        else
            c_i_Loop(i) = prodConc;
        end
        
    end
    
    LoopResults.ConcV     = LoopRes_ConcV;
    LoopResults.Prod_Conc = LoopRes_Prod_Conc;
    LoopResults.Pos_Conc  = LoopRes_Pos_Conc;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  C. CHECK IF ANY OF THE CONCENTRATIONS IS ABOVE OR     %
    %     BELOW THE LIMIT AND IF ALL REACTIONS ARE FEASIBLE  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Checking that all energetics are favourable
    [DGrVLoop]   = calcEnergetics(StM, Param, stoLoopM, n_ATPLoopV, protTranslocLoopV);
    DGCheck      = any(DGrVLoop > DGTol);
    
    %Check that DG has not changed between calculations
    sum_DG_check = sum_DG_loop - sum(DGrVLoop);
    sum_DG_loop  = sum(DGrVLoop);
    
    %Checking that all concentrations are withing physiological range
    concCheck    = or(any(c_i_Loop < Min_Conc), any(c_i_Loop > Max_Conc));
    

    %Check if first solution gives a feasible result (ignore C1 and Cn)
    posMinConcLoop = find(c_i_Loop(2:end-1) < Min_Conc,1) + 1;  %To compensate for ignoring C1
    if any(posMinConcLoop) 
        solutionFlag = -1;
        StM(LoopResults.Pos_Conc(posMinConcLoop)+1: id_Cout-1) = 1;
        LoopResults.Prod_Conc(posMinConcLoop+1 : end-1) = 1;
        %       LoopResults.Prod_Conc.(char(reacPath))(end) = Cn;
%         [DGrVLoop]   = calcEnergetics(StM, Param, stoLoopM, n_ATPLoopV, protTranslocLoopV);
        %       sprintf('No solution has been found in this configuration')
        break
    end
end
    
