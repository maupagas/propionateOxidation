%Define Table to plot results
refDG     = [-50, -55];    %kJ/mol
refRatioH = 9/3;    %Ratio H+/ATP
refCoA    = 0.001;  %CoA-SH
refpHin   = 7;      %pH
refPh2    = 1e-9;    %kJ/mol
refCO2    = 0.01;    %Ratio H+/ATP
refT      = 308.15;  %CoA-SH
refpHout  = 7;      %pH

%Create a combination matrix of all the results to plot
refPar = allcomb(refDG, refRatioH, refCoA, refpHin, refPh2, refCO2, refT, refpHout);

%Preallocate matrices to allow appendage of new results
totalfilteredVals = [];
idFilteredVals = [];

%For each row of the matrix to check,identify the values from the results matrix 
for i = 1:length(refPar(:,1))

%Find preselected values from the total inputs
filteredValues = pathwayInputs(pathwayInputs(:,1) == refPar(i,1) & pathwayInputs(:,2) == refPar(i,2) & pathwayInputs(:,3) == refPar(i,3) & ... 
                               pathwayInputs(:,4) == refPar(i,4) & pathwayInputs(:,5) == refPar(i,5) & pathwayInputs(:,6) == refPar(i,6) & ...
                               pathwayInputs(:,7) == refPar(i,7) & pathwayInputs(:,8) == refPar(i,8), :);

%Add previous results to new ones
totalfilteredVals = [totalfilteredVals; filteredValues];

%Find position of filteredValues
previdFiltVals = find(sum(refPar(i,:) == pathwayInputs,2) == length(pathwayInputs(1,:)));

%Add previous results to new ones
idFilteredVals = [previdFiltVals; idFilteredVals];
                                     
end
             

          