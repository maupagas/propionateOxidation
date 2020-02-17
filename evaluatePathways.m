%%%%%%%%%%%%%%%%% NEEDS FURTHER CLARIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a script written to calculate the concentration of the products
% in a given pathway. Concentrations of components are to be maintained
% within a range.

%It provides a different combination of proton translocation over the
%pathways in those reactions that are membrane bound.

%At the end of the script, reactions are evaluated, checking if the
%concentrations of the components are within the limits described, and all
%the reactions have a DG <0.

%% DEFINE PATHWAYS TO EVALUATE
% 
%Number of reactions to evaluate (To evaluate so far only the P1's (Reactions that do not contain any loop)
% firstPath  = 'P1a';
% lastPath   = 'P1b';
firstPath  = Reac.Pathway(1);
lastPath   = Reac.Pathway(end);
Param.firstPath2Eval = find(strcmp(char(firstPath), Reac.Pathway));
Param.lastPath2Eval  = find(strcmp(char(lastPath),  Reac.Pathway));
loopFlag = 0;       % Flag not to allow concentrations to go over the maximum limit (Only used in loop calculations)

%pKa's (to be adjusted with temperature) 
pKa_CO2 = 6.31; pKa_HCO3 = 10.26;

%% DEFINE DIFFERENT SCENARIOS TO EVALUATE THE PATHWAYS %%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% CREATE A NEW FUNCTION TO LOAD THIS #MPG: 12/9/19 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select the type of analysis to perform: 

%(1) Sensitivity 
%(2) All Possible Combinations
%(3) Hydrogen concentrations range to calculate methanogenic niche
combNames   = Param.combNames;
combinationAnalysis = char(Param.sensAnalysis);

switch combinationAnalysis
    case 'sensitivity'
        %Calculate Matrix Based on the Single Point Sensitivity Analysis
        [combValues] = calcPivotSensMatrix(Combination, Param.refComb, combNames);
        
    case 'allCombinations'
        % %Code for evaluate all combinations possible (Postprocessing of the data needs still to be managed)       
        combValues = allcomb(Combination.(char(combNames(1))), Combination.(char(combNames(2))), Combination.(char(combNames(3))), Combination.(char(combNames(4))), ...
                             Combination.(char(combNames(5))), Combination.(char(combNames(6))), Combination.(char(combNames(7))), Combination.(char(combNames(8))));
     
    case 'H2Range'

        %Define range of hydrogen concentrations
        H2_Range = [1.00E-10
                    1.00E-09
                    1.28E-08];
                
        %Overwrite Hydrogen Concentrations
        Combination.H2           = H2_Range;				

        %Calculate the "Sensitivity"Matrix
        [combValues] = calcPivotSensMatrix(Combination, Param.refComb, combNames);
        
end

        %DEBUGGING MODE (If it is >0, run the number of defined
        %combinations. If -1, define the custom condition to be evaluated
        numComb = 0;
        
        if numComb > 0
            combValues = combValues(1:numComb,:);
        elseif numComb == -1
            combValues = [-50	3.333333333  0.001	7	3.5E-9	0.01	308.15	7.465];
        end
        
        St.combinationList      = combValues;
        St.combinationListNames = combNames;


%Store list of Combination names
for k = 1:length(St.combinationList(:,1)), St.combinationNames{k} = sprintf('%s_%d','Combination',k); end
StM_0 = St.StM;
combinationListNames = St.combinationNames;

%Identifier for each one of the variables for the sensitivity analysis
id_DGATP = strcmp('DG_ATP', Param.combNames);
id_Ratio = strcmp('Ratio_H_ATP', Param.combNames);
id_CoASH = strcmp('CoA_SH', Param.combNames);
id_pHin  = strcmp('pH_in', Param.combNames);

id_H2    = strcmp('H2',     Param.combNames);
id_CO2   = strcmp('CO2',    Param.combNames);
id_T     = strcmp('T',      Param.combNames);
id_pHout = strcmp('pH_out', Param.combNames);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% UNTIL HERE THIS SHOULD GO IN ANOTHER FUNCTION %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Evaluate every possibility of DG_ATP, ratio_H_ATP, concCoA, concCO2
%% THIS IS THE MAIN FOR PATHWAY
for m = 1:length(combValues(:,1)) 
    
    sprintf('Calculating combination %d out of %d...', m, length(combValues(:,1)))
    
    %*************************************************************************
    %UPDATE COMBINATION VALUES: THIS HAS TO BE DEFINED BY HAND
    Param.DG_ATP = combValues(m, id_DGATP);          concH2  = combValues(m, id_H2);      
    ratio_H_ATP  = combValues(m, id_Ratio);          concCO2 = combValues(m, id_CO2); 
    concCoA      = combValues(m, id_CoASH);          T       = combValues(m, id_T); 
    pH_in        = combValues(m, id_pHin);           pH_out  = combValues(m, id_pHout);
    
    %Intracellular parameters
    Param.DG_Prot        = Param.DG_ATP / ratio_H_ATP;
    Param.n_H_ATP        = ratio_H_ATP;   
    St.StM(St.id.CoA_SH) = concCoA;
    St.StM(St.id.H)      = 10^-pH_in;
    
    %Extracellular parameters
    St.StM(St.id.H2)   = concH2;
    St.StM(St.id.CO2)  = concCO2;
    Param.T    = T; 
    St.StM(St.id.Hout) = 10^-pH_out;

    %**************************************************************************

    %UPDATE STATES IF REQUIRED
    %Update DG0f if T is increased - For those where Enthalpy is not
    %available, DG0ft is maintained as per DG0f
    G0ftM      = (Param.H0fM + (Param.T / 298.15)*(Param.G0fM - Param.H0fM));
    %Ignore NaNs (NEED TO READJUST STILL FOR CoA-bound components)
    G0ftM (isnan(G0ftM)) = Param.G0fM(isnan(G0ftM));
    Param.DG0ft = G0ftM;
    
    %Updates potential of the cell -> REMOVE THIS FUNCTION AS IT IS NOT
    %REQUIRED
    [Param.dp, Param.EV] = calc_dp(St, Param);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2.- Load electron Carrier concentrations     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % %Load Concentrations for electron Carriers as a function of the
    % %translocated protons (Assumed in equilibrium with H2)
    Reac.eC_Conc.ProtTransloc = -2:2;

    %Calculate electron carrier concentrations and trim it for the feasible
    %reactions
    [ratio_eCM, eC_RatioNames, protTranslocFeasComb, Res_eC_Conc, eCFeas_Conc, eC_ConcTable] = calcFeas_eC(Param, St, Reac);

    for nPath = Param.firstPath2Eval:Param.lastPath2Eval
        reacPath = Reac.Pathway(nPath);
        Results.(char(combinationListNames(m))).(char(reacPath)).ProtTranslocComb = protTranslocFeasComb.(char(reacPath));
        Results.(char(combinationListNames(m))).(char(reacPath)).eC_ratios        = ratio_eCM.(char(reacPath));
        Results.(char(combinationListNames(m))).(char(reacPath)).eC_Conc          = eCFeas_Conc.(char(reacPath));
        Results.(char(combinationListNames(m))).(char(reacPath)).eC_ReoxNames     = eC_RatioNames.(char(reacPath));
        Results.(char(combinationListNames(m))).(char(reacPath)).eC_ConcTable     = eC_ConcTable.(char(reacPath));
    end
    %Calculate DGr from excel for future calculations
    [DGr_FullSto] = calcEnergetics(St.StM, Param, Reac.stoM_xlsx,  0, 0);
      
    %*******************************************************************************

    %Evaluates each one of the Pathway (THIS IS THE MAIN FOR LOOP FOR EACH
    %PATHWAY)
    for k = Param.firstPath2Eval:Param.lastPath2Eval
        reacPath = Reac.Pathway(k);
        id_eC = find(strcmp(Reac.(char(reacPath)).Labels_ord, 'eC'), 1);
        id_eB = find(strcmp(Reac.(char(reacPath)).Labels_ord, 'eB'), 1);
        pathSolvingMode = Reac.(char(reacPath)).PathSolvingMode;
        
        sprintf('Evaluating Pathway %s ...', char(reacPath))
        
        %Load all elements required from the Reac structure
        stoM   = Reac.(char(reacPath)).stoM_ord;
        n_ATPV = Reac.(char(reacPath)).n_ATP_ord;
        
        %Matrices for constant and variable values, along with the position
        %of the concentration
        constConcM = Reac.(char(reacPath)).constConcM;
        varConcM   = Reac.(char(reacPath)).varConcM;
        posConcV   = Reac.(char(reacPath)).posConcV;
        
        %STORE ALL INPUT RESULTS: INPUTS, DGr, dp, DGProt
        Results.(char(combinationListNames(m))).(char(reacPath)).Inputs      = repmat(combValues(m,:), length(protTranslocFeasComb.(char(reacPath))), 1);
        Results.(char(combinationListNames(m))).(char(reacPath)).dp          = repmat(Param.dp,        length(protTranslocFeasComb.(char(reacPath))), 1);
        Results.(char(combinationListNames(m))).(char(reacPath)).DG_Prot     = repmat(Param.DG_Prot,   length(protTranslocFeasComb.(char(reacPath))), 1);
        Results.(char(combinationListNames(m))).(char(reacPath)).DGr_FullSto = repmat(DGr_FullSto,     length(protTranslocFeasComb.(char(reacPath))), 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        %Preallocate sizes
        ConcResM   = zeros(length(St.StM), length(Reac.(char(reacPath)).Names_ord));
        Prod_ConcM = zeros(length(protTranslocFeasComb.(char(reacPath))), length(Reac.(char(reacPath)).Names_ord));
        Pos_ConcV  = zeros(1, length(Reac.(char(reacPath)).Names_ord)); 
        DGrM       = Prod_ConcM;
        OptimConc  = zeros(length(St.StM), length(protTranslocFeasComb.(char(reacPath))));
        isProtTranslocFeasComb = isempty(protTranslocFeasComb.(char(reacPath)));
        eC_Path = fieldnames(eCFeas_Conc.(char(reacPath)));
        numReac = Reac.(char(reacPath)).numReac;
        %Collect how many loops are available in the pathway
        numLoops = idLoop.(char(reacPath)).numLoops;
        
        loopNameInit = 'Loop_1';
        calcLoops = 0;
        
        %Move This part of the code before doing it for each line of the code
        if numLoops > 0
            startLoop = idLoop.(char(reacPath)).(char(loopNameInit)).Start;
            endLoop   = idLoop.(char(reacPath)).(char(loopNameInit)).End;
                                    %Pass up the information required for the loop
            LoopData = idLoop.(char(reacPath)).(char(loopNameInit));                      
        else
            startLoop = 0;
            endLoop   = 0;
        end

    
        
        % Evaluates the pathway with all possible proton translocation combinations allowed  (This is the most consuming part of this script  
        for j = 1:length(protTranslocFeasComb.(char(reacPath))(:,1))
            %Recalls the Protons to be translocated for calculating Proton
            %Concentrations
            if isProtTranslocFeasComb
                continue
            else
                protTranslocV = protTranslocFeasComb.(char(reacPath))(j,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                %Update e-Carriers concentrations (either tabulated or calculated) 
                %THINK IF THIS CAN BE IMPROVED SOMEHOW (## MPG 28/11 ##)
                for i = 1:length(eC_Path)
                    eC = eC_Path(i);
%                     eC_reox = Param.Reox_reac(i);
%                     id_eC_reox = Reac.(char(reacPath)).id.(char(eC_reox)) == -1;
%                     if id_eC_reox
%                         continue
%                         % e-Carriers such as NADH and Fdred are directly reoxidized
%                         % to H2. Therefore, its concentrations are tabulated beforehand.      
%                     else
                        % e-Carriers that are not reoxidized directly to H2 need to be recalculated  
                        eCred_Conc = eCFeas_Conc.(char(reacPath)).(char(eC))(j);
                        St.StM(St.id.(char(eC))) =  eCred_Conc;
%                     end
                end
            end

            StM = St.StM;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% CALCULATION OF PRODUCTS CONCENTRATIONS

            % INITIALIZE FIRST %
            
            %Reactions counter and number of reactions in the while loop calculation
            counterRc = 1;
                
            %IF THERE IS A ELECTRON BIFURCATION, CALCULATE FIRST   
%             while pathSolvingMode == 2 && solutionFlag == 0  
            solutionFlag = 0;
            while solutionFlag ~= -1
                %If there is eB, solve first eB reaction
                if pathSolvingMode == 2
                    [ConcV, prodConc, r1, r2, eBreax] = solve_eB(Param, StM, St, Reac, stoM, n_ATPV, constConcM, varConcM, protTranslocV, reacPath);
                    
                    %Check if the ratio of the other eCarrier pair is within range
                    if r2 > 1e4 || r2 < 1e-4
                        break
                    else
                        StM = ConcV;
                    end
                end
                
                %Write the name of the loops for the first reaction
                if numLoops > 0
                    calcLoops = 0;
                    loopName = loopNameInit;
                else
                    loopName = [];
                end


                %Solve stoichiometry with loop
                while calcLoops <= numLoops && counterRc <= numReac

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %Move This part of the code before doing it for each line of the code
%                     if ~isempty(loopName) && calcLoops < numLoops
%                         startLoop = idLoop.(char(reacPath)).(char(loopName)).Start;
%                         endLoop   = idLoop.(char(reacPath)).(char(loopName)).End;
%                     else
%                         startLoop = 0;
%                         endLoop   = 0;
%                     end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %SOLVE PRODUCT SEQUENTIALLY
                    if counterRc < startLoop || counterRc > startLoop 

%                         sprintf('Solve sequentially reaction %d', counterRc)

                        %Update matrices
                        stoV  = stoM(:,counterRc);
                        n_ATP = n_ATPV(counterRc);
                        constConcV = constConcM(:,counterRc);
                        varConcV = varConcM(:,counterRc);
                        posConc = posConcV(counterRc);   %position of the product of the reaction

                        
                        % Solve sequentially: Calculate product concentration such as DGr = 0
                        [ConcV, prodConc] = calcProdConc(Param, StM, stoV, n_ATP, constConcV, varConcV, protTranslocV(counterRc), loopFlag);
                        
                        %Store Results into structures
                        ConcResM(:,counterRc)   = ConcV;
                        Prod_ConcM(j,counterRc) = prodConc;
                        Pos_ConcV(counterRc)    = posConc;
                        
                        %Updates the States Matrix with the calculated result
                        StM = ConcV;

                        counterRc = counterRc + 1;

                        %If any concentration is below the required concentration, the
                        %for loop is terminated
                        if prodConc < Param.Min_Conc && posConc ~= -1 
                            StM((posConc+1):(St.id.H2-1)) = 1;
                            solutionFlag = -1;
                            breakReac = counterRc;
                            break
                        end
                    
                    %SOLVE LOOP REACTION (either Ratio or loop block)
                    elseif counterRc == startLoop
                    % sprintf('Solve loop as a block reactions %d to %d', startLoop, endLoop)

                        
                        %Differentiate if it is a ratio calculation
                        %reaction or a Loop solving reaction
                        isRatioRc = strcmp(Reac.(char(reacPath)).Labels_ord(counterRc), 'Cr');
                        if isRatioRc
                            stoV  = stoM(:,counterRc);
                            n_ATP = n_ATPV(counterRc);
                            constConcV = constConcM(:,counterRc);
                            numProtTransloc = protTranslocV(counterRc);
                            ratioFlag = 1;  %Loop solving Ratio Flag 
                            [ConcV, prodConc, posConc, ratiospc] = solveRatio(StM, stoV, n_ATP, constConcV, numProtTransloc, Param, LoopData);
%                             sprintf('Calculated ratio of species is %d', ratiospc)
                            
                            %Store Results into structures
                            ConcResM(:,counterRc)   = ConcV;
                            Prod_ConcM(j,counterRc) = prodConc; 
                            Pos_ConcV(counterRc)    = posConc;
                            
                            %Update Loop counter                                
                            counterRc = counterRc + 1;
                            
                            %If the ratio is below or above the required
                            %value the while loop is terminated
                            if ratiospc < 1e-4 || ratiospc > 1e4 
                                StM((posConc+1):(St.id.H2-1)) = 1;
                                solutionFlag = -1;
                                breakReac = counterRc;
                                break
                            end
                            
                        else                 
                            %Solve Loop
                            protTranslocLoopV = protTranslocFeasComb.(char(reacPath))(j, startLoop:endLoop);
                            [ConcV, LoopResults, solutionFlag] = solveLoop(StM, Param, LoopData, protTranslocLoopV, j);
                            
                            %Store Results into structures
                            ConcResM(:, counterRc:endLoop)   = LoopResults.ConcV;
                            Prod_ConcM(j, counterRc:endLoop) = LoopResults.Prod_Conc; 
                            Pos_ConcV(:, counterRc:endLoop)  = LoopResults.Pos_Conc;
                            
                            if solutionFlag == -1
                                breakReac = counterRc;
                                StM = ConcV; 
                                break
                            end
                            
                            %Update Loop counter                                
                            counterRc = endLoop + 1;
                        end
                        
                        %Updates the States Matrix with the calculated result
                        StM = ConcV;                              

                        % Update number of loops calculated counter
                        calcLoops = calcLoops + 1;
                        if calcLoops < numLoops
                            loopName = sprintf('Loop_%d', calcLoops+1);    
                            startLoop = idLoop.(char(reacPath)).(char(loopName)).Start;
                            endLoop   = idLoop.(char(reacPath)).(char(loopName)).End;
                            LoopData = idLoop.(char(reacPath)).(char(loopName));
                        end
                    end
                end

                % sprintf('Whole stoichiometry with loop has been solved')       
                %Check Gibbs reactions
                if pathSolvingMode == 2 
                    [DGr] = calcEnergetics(StM, Param, stoM, n_ATPV, protTranslocV);
                    breakReac = find(DGr>0,1);

                    if all(DGr > 1e-4)
                        solutionFlag = -1;
                    elseif all(DGr < 1e-4)
                        solutionFlag = 1;
                    end
                end

                %Terminate the while loop for the electron bifurcation (in case
                %there is none)
                if solutionFlag == 0
                    if pathSolvingMode ~= 2
                        break
                    elseif pathSolvingMode == 2 && eBreax(breakReac) == 0
%                         sprintf(' The calculation of the eB reaction does not have a solution');
                        break
                    elseif pathSolvingMode == 2 && eBreax(breakReac) == 1

                        %Ferredoxin ratio is recalculated
                        stoV  = stoM(:,breakReac);
                        n_ATP = n_ATPV(breakReac);
                        constConcV = ones(size(stoV));
                        numProtTransloc = protTranslocV(breakReac);
                        constConcV(St.id.Fdred) = 0;

                        %Calculate new ferredoxin ratio
                        [ConcV, prodConc, posConc, ratiospc] = solveRatioFerredoxin(StM, stoV, n_ATP, constConcV, numProtTransloc, Param, St);
                        StM = ConcV;
                        solutionFlag = 0;
                        %If the ratio is below or above the required
                        %value the while loop is terminated
                        if ratiospc < 1e-4 || ratiospc > 1e4 
                            solutionFlag = -1;
%                             sprintf('Ferredoxin ratio is out of the boundaries and calculation is terminated')
                            break
                        end
                    end
                end
                %Restart reaction break counter
                breakReac = 0;
                if solutionFlag == 1, break, end
            end
            %Collects the concentrations optimized for minimizing dissipated DG
            %in each reaction and resets the state matrix for the next mode calculation

            OptimConc(:,j) = StM;   % NEED TO PREALLOCATE THIS VARIABLE

            %States Matrix is re-started for the next calculation 
            StM = St.StM;

            %Calculate Energetics of the reactions 
            a_i    = OptimConc(:,j);
%             reacM  = Reac.(char(reacPath)).stoM_ord;
            n_ATPV = Reac.(char(reacPath)).n_ATP_ord;

           [DGrV] = calcEnergetics(a_i, Param, stoM, n_ATPV, protTranslocV);

           %Find the first positive DG value (0.01 used to ignore values that are
           %bigger than 0 but negligible -e.g. 1e-14)
%                posDG_positive = [ find(DGrV > 0.01, 1), find(Prod_ConcM(j,:) < Param.Min_Conc, 1)] ;
           posDG_positive = find(DGrV > 0.01, 1) ;

           posDG_positive = min(posDG_positive) + 1; 

           %Correct those values with positive Gibbs to values of 1000
           if numel(posDG_positive)>0
               if isempty(id_eB)
                   DGrV((posDG_positive(1)):(id_eC -1)) = 1000;
               else
                   DGrV((posDG_positive(1)):(id_eB -1)) = 1000;
               end
           end
           
           DGrM(j,:) = DGrV;
           
        end

           %Store all results in Structure
            Results.(char(combinationListNames(m))).(char(reacPath)).ConcV     = ConcResM;
            Results.(char(combinationListNames(m))).(char(reacPath)).Prod_Conc = Prod_ConcM; 
            Results.(char(combinationListNames(m))).(char(reacPath)).Pos_Conc  = Pos_ConcV;
            Results.(char(combinationListNames(m))).(char(reacPath)).DGrV      = DGrM;
            Results.(char(combinationListNames(m))).(char(reacPath)).OptimConc = OptimConc;   % NEED TO PREALLOCATE THIS VARIABLE


    end

    %% Check both Concentrations and DG Compatibilities
    for i = Param.firstPath2Eval:Param.lastPath2Eval

        reacPath = Reac.Pathway(i);
        if ~isempty(protTranslocFeasComb.(char(reacPath)))
            %If reaction is bigger or lower than Max concentration, is labeled as
            %1, if it is not, is labeled as 0
            Results.(char(combinationListNames(m))).(char(reacPath)).ConcCheckM = and(or(Results.(char(combinationListNames(m))).(char(reacPath)).OptimConc(1:St.id.H2-1,:) > Param.Max_Conc, ...
                                                                                         Results.(char(combinationListNames(m))).(char(reacPath)).OptimConc(1:St.id.H2-1,:) < Param.Min_Conc), ...
                                                                                         Results.(char(combinationListNames(m))).(char(reacPath)).OptimConc(1:St.id.H2-1,:)    ~= 1);

            Results.(char(combinationListNames(m))).(char(reacPath)).DGCheckM = Results.(char(combinationListNames(m))).(char(reacPath)).DGrV <= 1e-6;

            %For having all concentration within the required range, all the values
            %for each COLUMN in the checking matrix have to be 0 ( BETTER TO USE AN ANY)
            Results.(char(combinationListNames(m))).(char(reacPath)).Compatible_Conc = all(Results.(char(combinationListNames(m))).(char(reacPath)).ConcCheckM == 0);
            %For having all reactions under a negative DG, all the values for each ROW in the checking matrix have to be 1
            Results.(char(combinationListNames(m))).(char(reacPath)).Compatible_DG   = all(Results.(char(combinationListNames(m))).(char(reacPath)).DGCheckM' == 1);
            %Checks reactions that are compatible in both Concentrations and Gibbs Free Energy
            Results.(char(combinationListNames(m))).(char(reacPath)).Compatible_Combination = and(Results.(char(combinationListNames(m))).(char(reacPath)).Compatible_Conc == 1, Results.(char(combinationListNames(m))).(char(reacPath)).Compatible_DG == 1);
            %Finds the reaction that satisfies both compatibilities and returns a message
            Results.(char(combinationListNames(m))).(char(reacPath)).CompCombId = find(Results.(char(combinationListNames(m))).(char(reacPath)).Compatible_Combination == 1);

            %Group results for group of net H translocation
       %     numNetHTransloc = Reac.ProtTranslocComb.(char(reacPath))(Results.(char(combinationListNames(m))).CompCombId.(char(reacPath))) + Param.num_H_ATP;

            netHTranslocation    = sum(protTranslocFeasComb.(char(reacPath)), 2) + Param.n_ATP_SLP * Param.n_H_ATP;
            netATP_Prod = netHTranslocation / Param.n_H_ATP;
            netDG_Prod  = netHTranslocation * Param.DG_Prot;
            numValidNetHTransloc = netHTranslocation(Results.(char(combinationListNames(m))).(char(reacPath)).CompCombId);

            
            for w = min(numValidNetHTransloc):max(numValidNetHTransloc)

                %Name the character
        %         HTranslocName = sprintf('NetHTransloc_%s',num2str(w));
        %         HTranslocName = strrep(HTranslocName,'-','_');
                id_NetH = find(Param.netHTranslocRange == round(w,2));
                netHTranslocName = Param.netHTranslocNames(id_NetH);

                %Valid results combination
                validResults = Results.(char(combinationListNames(m))).(char(reacPath)).CompCombId;
                %Group results for each netproton translocation
                Results.(char(netHTranslocName)).(char(reacPath)).(char(combinationListNames(m))) =  validResults(round(numValidNetHTransloc,2) == round(w,2));

            end
            
        else        
            netHTranslocation    = [];
            netATP_Prod = [];
            netDG_Prod  = [];
           
        end
        
            Results.(char(combinationListNames(m))).(char(reacPath)).netHTransloc = netHTranslocation;
            Results.(char(combinationListNames(m))).(char(reacPath)).netATP_Prod  = netATP_Prod;
            Results.(char(combinationListNames(m))).(char(reacPath)).netDG_Prod   = netDG_Prod;
    %     
    %     if sum(Compatible_Conc > 0),
    %         sprintf('The concentrations combination of %d is compatible with the cell requirements for Mode %s', Compatible_Conc, char(Reac.Pathway(i)))
    %     else 
    %         sprintf('NO COMBINATION OF CONCENTRATIONS have been found compatible in this Metabolic Route %s', char(Reac.Pathway(i)))
    %     end
    end

end

%Collect ONLY the results for the combinations that are feasible
[PrintResults] = printFeasibleResults(Results, Reac, Param, combinationListNames);
  
clearvars -except Reac Param idLoop St Results PrintResults Combination