% This function builds a matrix with single point sensitivity analysis from
% the parameters defined in the spreadsheet. At the end of the matrix
% built, it deletes those ones that are repeated.

function [combinationMatrixFinal] = calcPivotSensMatrix(Combination, refComb, combNames)

%Calculate the possible number of single-point sensitivity analysis
%combinations
for i = 1:length(combNames)
    numComb(i) = numel(Combination.(char(combNames(i))));
end

%Sum all the number of possible combinations
numPossComb = sum(numComb);

%Create the combination matrix for all possibilities
combinationMatrix = zeros(numPossComb, length(refComb));

%Collect the number of accumulated variables for each parameter to be
%studied, in order to build the sensitivity analysis matrix
numCumComb = cumsum(numComb);

%Initialize the position of the value and the value of the variable used in
%the previous iteration
idPosVal = 0;     idVarPrev = 0;

%Fill the matrix with all the possible combinations
for i = 1:numPossComb
    
    %Identify to which variable belongs the combination 
    %(e.g. iterations 1-6 for DG_ATP, 7-13 to Ratio H+/ATP, etc....
    idVar = find(i <= numCumComb,1);
    
    %Checks if the variable has changed or not
    diffIdVar = idVar - idVarPrev;
    %Store the column of the variable for the next iteration
    idVarPrev = idVar;
    
    %Find position between vectors: If we are using the same variable, 
    %pass onto the next value. If there variable changes, re-start at 1.
    if diffIdVar == 0
        idPosVal = idPosVal + 1;
    else 
        idPosVal = 1;
    end
    
    %Fill single-point sensitivity matrix:
    % 1) Fill the row with the reference vector 
    combinationMatrix(i, :)     = refComb;
    % 2) Change the variable of the reference vector with the new variable
    combinationMatrix(i, idVar) = Combination.(char(combNames(idVar)))(idPosVal);

end

%Delete all the rows that are repeated (e.g. Reference value rows)
combinationMatrixFinal = unique(combinationMatrix, 'rows');

