%
combination = 'H2Range';

function [Combination] = loadCombsToEval(combination, 

switch combination
    case 'sensitivity'
        Combination.varNames = {'DG_ATP', 'ratio_H_ATP', 'CoA_SH', 'pH_in', 'H2', 'CO2', 'T', 'pH_out'}; 
        
        Combination.DG_ATP       = [-40	        -45         -50	        -55	        -60	      -65];	
        Combination.ratio_H_ATP  = 9/3:1/3:15/3;
        Combination.CoA_SH       = [1.00E-06	1.00E-5     1.00E-04	1.00E-03	1.00E-02];			
        Combination.pH_in        = [6           6.5         7           7.5         8];				
        Combination.H2           = [1.00E-10    1.00E-09	1.28E-08	2.50E-07];				
        Combination.CO2	         = [1.00E-06	1.00E-05    1.00E-04	1.00E-03    0.01];				
        Combination.T	         = [25	        35	        45          55] + 273.15;				
        Combination.pH_out       = [6	        7	        8];		
        
        %Provide reference vector for each variable (KEEP IN ORDER)
%         refV = [4, 2, 4, 3, 2, 5, 2, 2];
        refV = [3,  ... %DG_ATP 
                1,  ... %Ratio H/ATP
                4,  ... %[CoA-SH] 
                3,  ... %pH_in
                2,  ... %[H2] 
                5,  ... %[CO2] 
                2,  ... %T 
                2]; ... %pH_out 
                        
        [combValues, refComb] = calcPivotSensMatrix(Combination, refV);
        
         Param.refComb = refComb;
%         combinationListNames = Combination.varNames;        
        
    case 'all'
        % %Code for evaluate all combinations possible (NOT TO BE USED YET)
        DG_ATP  = [-65, -60];
        ratio_H_ATP = [9/3 10/3 11/3 12/3];
        concCoA = 0.01;
        pH_in   = 7;
        % ***** EXTRACELLULAR
        conc_H2  = 1.29e-8;
        conc_CO2 = 1e-4;
        T        = 308.15;
        pH_out   = 7;
           
        combValues      = allcomb(DG_ATP, ratio_H_ATP, concCoA, pH_in, conc_H2, conc_CO2, T, pH_out);
        
    case 'H2Range'
        Combination.varNames = {'DG_ATP', 'ratio_H_ATP', 'CoA_SH', 'pH_in', 'H2', 'CO2', 'T', 'pH_out'}; 

        H2_Range = [1.00E-10
                    1.00E-09
                    1.28E-08];

        
        Combination.DG_ATP       = -50;	
        Combination.ratio_H_ATP  = 10/3;
        Combination.CoA_SH       = 1.00E-03;			
        Combination.pH_in        = 7;				
        Combination.H2           = H2_Range;				
        Combination.CO2	         = 0.01;				
        Combination.T	         = 35 + 273.15;				
        Combination.pH_out       = 7;   %pH BSM2;	
        
        refV = [1, ... %DG_ATP 
        1,  ... %Ratio H/ATP
        1,  ... %[CoA-SH] 
        1,  ... %pH_in
        2,  ... %[H2] 
        1,  ... %[CO2] 
        1,  ... %T 
        1]; ... %pH_out 
                        
        [combValues, refComb] = calcPivotSensMatrix(Combination, refV);
        
        Param.refComb = refComb;
        
end
        combValuesNames = {'DG_ATP', 'ratio_H_ATP', 'CoA_SH', 'pH_in', 'H2', 'CO2', 'T', 'pH_out'};

        %DEBUGGING MODE (If it is >0, run the number of defined
        %combinations. If -1, define the custom condition to be evaluated
        numComb = 0;
        
        if numComb > 0
            combValues = combValues(1:numComb,:);
        elseif numComb == -1
            combValues = [-50	3.333333333  0.001	7	3.5E-9	0.01	308.15	7.465];
        end
        
        St.combinationList      = combValues;
        St.combinationListNames = combValuesNames;


%Store list of Combination names
for k = 1:length(St.combinationList(:,1)), St.combinationNames{k} = sprintf('%s_%d','Combination',k); end
StM_0 = St.StM;
combinationListNames = St.combinationNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% UNTIL HERE THIS SHOULD GO IN ANOTHER FUNCTION %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%