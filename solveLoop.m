% A loop is defined by the following reaction
% Cin + Cn  <----> Cout + C1
% DGr = 0, Cout = 1e-6 M for calculating C1max and minimum ratio

function [StM, LoopResults, solutionFlag] = solveLoop(StM, Param, LoopData, protTranslocLoopV, nProtComb)

%Preload States
DG0ft    = Param.DG0ft;        Rth      = Param.Rth;          T  = Param.T;        
DG_Prot  = Param.DG_Prot;      DG_ATP   = Param.DG_ATP;
Max_Conc = Param.Max_Conc;     Min_Conc = Param.Min_Conc;

%Load matrices
stoLoopM   = LoopData.stoLoopM;      
n_ATPLoopV = LoopData.n_ATPLoopV;    
constConcLoopM = LoopData.constConcLoopM;
varConcLoopM   = LoopData.varConcLoopM;
posConcLoopV       = LoopData.posConcLoopV;
% protTranslocLoopV = LoopData.protTranslocLoopV(nProtComb,:);

%Identifiers of elements of the loop
id_C1 = LoopData.id_C1;              id_Cout = LoopData.id_Cout;
id_Cn = LoopData.id_Cn;
% startLoop = LoopData.StartLoop;
%
nonC1V   = [1:id_C1-1   id_C1+1:length(DG0ft)];
nonCoutV = [1:id_Cout-1 id_Cout+1:length(DG0ft)];


%INITIALIZE FIRST WHILE LOOP
concCheck = 1;         DGCheck = 1;   ratioCheck = 0;
ratioNew = 0;          ratioFlag = 0;
sum_DG_loop = 1000;    sum_DG_check = 1000;
Cn = Max_Conc;         Cout = Min_Conc;
DGTol = 1e-6;
ratiobiggerFlag = 0;

%Concentration of Cout is the lowest possible to give 'all the energy
%available' to the Loop
StM(id_Cout) = Cout;

%Preallocate LoopResults
LoopRes_ConcV     = zeros(length(StM), length(stoLoopM(1,:)));
LoopRes_Prod_Conc = zeros(1,length(stoLoopM(1,:)));
LoopRes_Pos_Conc  = LoopRes_Prod_Conc;

%% Calculate concentration of first component
while concCheck == 1 || DGCheck == 1 || ratioCheck == 0
    %         if ratioNew == 0,  sprintf('First calculation is being performed. \n'), else sprintf('Ratio is smaller than expected ratio. \n', sum(DGrVLoop)), end
    c_i_Loop = 0 * stoLoopM(1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%  CALCULATE LOOP CONCENTRATIONS  %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:length(stoLoopM(1,:))
        stoV = stoLoopM(:,i);
        ObjDG = 0 - protTranslocLoopV(i) * Param.DG_Prot;
        n_ATP_Loop = n_ATPLoopV(i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % A. SOLVE FIRST REACTION BY PROVIDING A VALUE FOR C1 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i == 1
            
            %If it is the first calculation, no new ratio is provided
            if ratioNew == 0
                
%                 %Calculates the Product concentration for a DG = 0
%                 prodStoich = sum(stoV .* varConcV);
%                 prodConc = exp(((-sum_u_i/prodStoich - DG0ft_var) / (Param.Rth * Param.T )));
                
                %Calculate product assuming a ratio (PREALLOCATE id_C1
                u_i = DG0ft(nonC1V) + Rth * T * log(StM(nonC1V));
                sum_u_i = sum(u_i .* stoV(nonC1V)) + ObjDG - DG_ATP * n_ATP_Loop;
                
                %Calculates the Product concentration (and the initial ratio) for a DG = 0
                C1 = exp((-sum_u_i/stoV(id_C1) - DG0ft(id_C1)) / (Rth * T));
                
                %Define the minimum ratio (if no energy is dissipated)
                ratio_min = Cn/C1;
                
                %If C1 is bigger than the maximum allowed concentration,
                %the reaction loop has to dissipate energy (If DG = 0 with
                %C1 > Max_Conc, it will be also favourable if C1 <
                %Max_Conc, as C1 is a product of the loop reaction
                if C1 > Max_Conc 
                    C1 = Max_Conc; 
                elseif C1 < Min_Conc
                    break
                end
          
            else    %For the rest of calculations in the loop, C1 concentration is provided
                C1 = C1_new;
            end
            
            %Update states vector
            StM(id_C1) = C1;
            prodConc = C1;
            LoopRes_Prod_Conc(i) = prodConc;
            
            %Update new ratio
            if ratioNew ~= 0
                ratio = ratioNew;
            end
            %             sprintf('Ratio of concentrations in the loop is %.2d. \n', ratio);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % B. CALCULATE SEQUENTIALLY REST OF THE PRODUCT CONCENTRATIONS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            %To allow the concentration of Cn to go over the physiological
            %limits only when the ratio has not been determined that it
            %will be bigger
            if i == length(stoLoopM(1,:)) && ratioFlag < 2
                loopFlag = 1;
            else
                loopFlag = 0;
            end
            
            %Calculate product concentration of the rest of the components
            
            constConcV = constConcLoopM(:,i);
            varConcV = varConcLoopM(:,i);
            posConc = posConcLoopV(i);

            [ConcV, prodConc] = calcProdConc(Param, StM, stoV, n_ATP_Loop, constConcV,  varConcV, protTranslocLoopV(i), loopFlag);
                        
            %Store Results into structures
            LoopRes_ConcV(:,i)   = ConcV;
            LoopRes_Prod_Conc(i) = prodConc;
            LoopRes_Pos_Conc(i)  = posConc;
            
            %If product concentration of Cn is below the concentration limit,
            %DO NOT UPDATE states matrix.
            if i == length(stoLoopM(1,:)) && prodConc < Min_Conc
                continue
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
    
%     if any(c_i_loop < Min_Conc) && ratioFlag ~= 2
%         solutionFlag = -1;
%         sprintf('No solution has been found in this configuration')
%     end
    %Check if first solution gives a feasible result (ignore C1 and Cn)
    posMinConcLoop = find(c_i_Loop(2:end-1) < Min_Conc,1) + 1;  %To compensate for ignoring C1
    if any(posMinConcLoop) && ratioFlag ~= 2
        solutionFlag = -1;
        StM(LoopResults.Pos_Conc(posMinConcLoop)+1: id_Cn-1) = 1;
        StM(id_Cn) = Cn;
        LoopResults.Prod_Conc(posMinConcLoop+1 : end-1) = 1;
        %       LoopResults.Prod_Conc.(char(reacPath))(end) = Cn;
        [DGrVLoop]   = calcEnergetics(StM, Param, stoLoopM, n_ATPLoopV, protTranslocLoopV);
        %       sprintf('No solution has been found in this configuration')
        break
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % D. PROVIDE A NEW VALUE FOR C1 TO RUN THE LOOP AGAIN (IF NECESSARY)  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calculate new ratio
    Cn_new   = prodConc;
    ratioNew = Cn_new / C1;
    
    %Compare the new ratio with the originla ratio
    if ratioNew < ratio_min
        ratioComparison ='smaller';
        ratioFlag = 1;
    elseif ratioNew >= ratio_min
        ratioComparison = 'bigger';
        ratioFlag = 2;
    end
    
    %C1 calculation is approach in different manners as function of the
    %ratio being bigger or smaller.
    switch ratioComparison
        case 'smaller'
            
            %A. Check that a solution has been found (regarding DG)
            if DGCheck == 1
                solutionFlag = -1;
                %                  sprintf('Ratio is smaller and a reaction with positive DG has been found. Therefore, loop calculation is terminated')
                break,
            elseif DGCheck == 0 && sum_DG_check == 0
                ratioCheck = 1;  %A solution has been found
            else
                %B. If all reactions had DG negative, a new concentration might be found
                C1_new = Cn_new / ratio_min;
                if C1_new < Min_Conc
                    StM(id_C1) = C1_new;
                    StM(id_Cn) = Cn_new;
                    solutionFlag = -1;
                    %                 sprintf('Loop calculation terminated. \n. Ratio is smaller than the minimum and product concentration falls below physiological concentrations' \n)
                    break
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CALCULATE THE CASE FOR A BIGGER RATIO  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
        case 'bigger'
         
            %IfCn is bigger than the physiological concentrations. Therefore,
            %we need to put C1 back into feasible concentrations.
            if Cn_new > Max_Conc, Cn_new = Max_Conc; end 
            StM(id_Cn) = Cn_new;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  A. Ratio is bigger and sum of DG is 0  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            if sum_DG_loop >= -DGTol && sum_DG_loop <= DGTol

                C1_new = Cn_new / ratioNew;
                %Check in the second run of the loop
                if ratiobiggerFlag == 1, ratioCheck = 1; end
                ratiobiggerFlag = 1;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  B. Ratio is bigger and sum of DG is negative %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            elseif sum_DG_loop <= -DGTol
                
                %Check in the second run of the loop
                if  ratiobiggerFlag == 2
                    ratioCheck = 1; 
                else
                    ratiobiggerFlag = 2;
                end
                %Allocates ALL THE ENERGY from the reactions in the loop to
                %the first loop reaction (2 reactants and two products)
                u_i_loop = DG0ft(nonC1V) + Rth * T * log(StM(nonC1V));

                %Sums all the chemical potentials except the species that is going to
                %change. Also adds the equivalent Gibbs of the translocated protons
                sum_u_i_loop = sum(u_i_loop .* stoLoopM(nonC1V,1)) - sum_DG_loop - DG_ATP * n_ATPLoopV(1);    
                
                %Calculated the new theoretical C1 for the loop reaction
                C1_new = exp(((-sum_u_i_loop - DG0ft(id_C1)) / (Rth * T)));

                if C1_new < Min_Conc
                    ratioNew = ratioNew * C1_new/Min_Conc;
                    C1_new = Cn_new /ratioNew;
                end
                
                if C1_new > Max_Conc 
                    C1_new = Max_Conc;             
                end              
                
            end
            
            %If any concentration falls below the minimum, ratio needs to
            %be recalculated
            if any(c_i_Loop < Min_Conc)
                ratioNew = ratioNew * min(c_i_Loop)/Min_Conc;
                C1_new = Cn_new / ratioNew;
                if C1_new > Max_Conc 
                    C1_new = Max_Conc;             
                end         
            end
            StM(id_C1) = C1_new;
%             StM(id_Cn) = Cn_new;
       
    end
    
    %A solution has been found
    solutionFlag = 1;
    
end

%If a solution has been found, calculate Cout
if solutionFlag == 1
    %%ratioFinal = ratioNew;
    %Calculate Cout
    stoV = stoLoopM(:,1);
    
    %Calculate product assuming a ratio
    ObjDG = 0 - protTranslocLoopV(1) * DG_Prot;
    
%     u_i = DG0f(nonCoutV) + Rth * T * log(StM(nonCoutV)) + chrV(nonCoutV).* F  .* EV(nonCoutV);
    u_i = DG0ft(nonCoutV) + Rth * T * log(StM(nonCoutV));
    
    sum_u_i = sum(u_i .* stoV(nonCoutV)) + ObjDG - DG_ATP * n_ATP_Loop(1);
    
    %Calculates the Product concentration for a DG = 0
%     Cout = exp(((-sum_u_i - DG0f(id_Cout) - chrV(id_Cout) * F .* EV(id_Cout)) / (Rth * T )));
    Cout = exp(((-sum_u_i/stoV(id_Cout) - DG0ft(id_Cout)) / (Rth * T )));

    if Cout > Max_Conc, Cout = Max_Conc; end
    StM(id_Cout) = Cout;
    
    %Check that only the first reaction changed its DG
    [DGrVLoop]   = calcEnergetics(StM, Param, stoLoopM, n_ATPLoopV, protTranslocLoopV);
    
    %Find the first positive DG value (0.01 used to ignore values that are bigger than 0 but negligible -e.g. 1e-14)
    posDG_positive = find(DGrVLoop > 1e-6, 1);
    
    if numel(posDG_positive)>0
        DGrVLoop(posDG_positive(1):end) = 1000;
    end
    
    
end


%%%% FOR DEBUGGING PURPOSES  %%%%%%%%%%
% if mod(nProtComb, 100) == 0
%     sprintf("Iteration %d has been solved", nProtComb)
% end