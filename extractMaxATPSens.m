%%%This script is used to extract the data from the sensitivity analysis
%%%for the different physiological and environmental conditions. 
%%% 
%%% The code still needs to be prepared to extract all possible
%%% combinations (21/11/19)

%Collect the variables and the reference values that need 
%to be extracted (the same as the sensitivity analysis) 
var2Extract = Param.combNames;
refPoint    = Param.refComb;

%List of pathways to assess
reacPathway = Reac.Pathway(Param.firstPath2Eval:Param.lastPath2Eval);

%Evaluate each variable
for m = 1:length(var2Extract)
    
    %Identify the variable from which extract the results
    varResults = var2Extract(m);
    %Preallocate the names for the different values used for each variable
    varValuesNames = cell(1,length(Combination.(char(varResults))));
    %Find the column of the variable to extract the values
    var_col = strcmp(var2Extract, varResults);
    %Obtain all values of the variable to extract its values
    varValues = Combination.(char(varResults));
    Output.varValues.(char(varResults)) = varValues;
    
    %Write names of the different number of Outputs (Think also if it is
    %better to give them a more significant name rather than 1, 2 3)
    for i = 1:length(varValues)
        varValuesNames{i} = sprintf(strcat(char(varResults), '_%d'), i);
    end
    
    %Preallocate matrix of all results for all variables for max ATP for
    %each combination
    plotResults = zeros(length(varValuesNames), length(reacPathway));

    
    %For each one of the variables, obtain the maximum yield of ATP for each pathway  
    for k = 1:length(varValuesNames)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    THIS EVENTUALLY SHOULD BE A SWITCH CASE IF IT IS A      %
        %    SENSITIVITY ANALYSIS OR AN ALL POSSIBILITIES SCENARIO   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Evaluate the sensitivity point corresponding to the variable evaluated
        sensPointV = refPoint .* (var_col==0) + varValues(k)* var_col;
                
        %Find the vector inside the huge sensitivity matrix 
        checkupVal = bsxfun(@minus,St.combinationList, sensPointV);      
        combId = find(all(checkupVal' == 0));        %if you deduct one from another should be all zeros
        combNum = sprintf('Combination_%d', combId); %Finds the combination number
        
        % For each pathway
        for j = 1:length(reacPathway)
            reacPath = reacPathway(j);
            
            %If there is any compatible combination, obtain the ATP produced and the Gibbs Free energy of the reaction
            if isfield(Results.(char(combNum)).(char(reacPath)), 'Compatible_Combination')
                feasibleComb = (Results.(char(combNum)).(char(reacPath)).Compatible_Combination);     %Compatible combinations are those with positive ATP and concentrations within the range
                feasATP_Prod =  Results.(char(combNum)).(char(reacPath)).netATP_Prod(feasibleComb);   %Collects the ATP of the feasible combinations
                Output.DGr.(char(varValuesNames(k))) = round(Results.(char(combNum)).(char(reacPath)).DGr_FullSto(1),0); %The Gibbs energy of the reaction is obtained for the plots

                %If ATP >= 0 exist, collect the maximum ATP obtained in the
                %given pathway for the given combination. If there are ATPs
                %lower than 0, puts a -10 as a dummy value 
                if  ~isempty(feasATP_Prod)
                    Output.Results.(char(varValuesNames(k))).max_ATP(j) =  max(feasATP_Prod);
                else
                    Output.Results.(char(varValuesNames(k))).max_ATP(j) =  -10;
                end 
                    
            %Define also as -10 if there are no compatible combinations    
            else
                Output.Results.(char(varValuesNames(k))).max_ATP(j) =  -10;
            end          
        end
   
        %Register results as matrices for bar plots        
        plotResults(k,:) = Output.Results.(char(varValuesNames(k))).max_ATP;
            
    end 
    
        %The vector of results for the variable is stored under the Output structure
        Output.plotResults.(char(varResults)) = plotResults;

end

clearvars -except Reac Param idLoop St Results PrintResults Combination Output
