% Identifies the substrates and the products for each reaction in each pathway 

function [constConcM, varConcM, posConcV, idCalcReax, numReac, solvePathwayID] = idProdsSubs(Reac, St, reacPath)

    stoM   = Reac.(char(reacPath)).stoM_ord;
    id_H2  = St.id.H2;

    %Identifies the components that participate in the reaction
    Conc2Change = (stoM ~= 0);
    Subs_id = (stoM == -1); %Identify substrates
    % Conc2Change = Conc2Change .* (Subs_id == 0);      % Keeps the concentration of the substrate also constant
    Conc2Change = bsxfun(@times,Conc2Change, (Subs_id == 0));      % Keeps the concentration of the substrate also constant
    Conc2Change(id_H2:end,:) = 0;                    % Keeps the concentrations of carriers, ATP, Pi, and CO2 as constant (for now)
    
    %Keep parameters labeled as constant as 0 in the matrix concentration
    %of variables changing
    idConstConc = strcmp(St.StType, 'c');
    Conc2Change(idConstConc,:) = 0;
    
    constConcM = (Conc2Change==0);
    varConcM   = (Conc2Change==1);
    [posConcV, ~] = find(varConcM > 0);
  
    %Identify Loops and electron bifurcation reactions
    id_eC = find(strcmp(Reac.(char(reacPath)).Labels_ord, 'eC'), 1);
    id_eB = find(strcmp(Reac.(char(reacPath)).Labels_ord, 'eB'), 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Identify if the reaction contains a loop or not %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Identify substrates and products for each step in the pathway
    id_Metabolites = (stoM ~= 0); 
    id_Metabolites(St.id.H2:end,:) = 0;             % To consider that the values of eC, ATPs and other metabolites (i.e. CO2) remain constant
    numSubs  = sum(id_Metabolites .* stoM < 0);
    numProds = sum(id_Metabolites .* stoM > 0);

    %Identify reactions with more than one substrate and one product
    idCalcReax = find(and(numSubs > 1, numProds > 1));

    %Find the number of reactions to be evaluated
    numReac = find(strcmp(Reac.(char(reacPath)).Labels_ord, 'Fp')) - 1;

    %Identify solving mode (1: Linear, 2: Loop solving, 3: Loop and
    %electron bifurcation
    if isempty(idCalcReax)
        solvePathwayID = 0;
    elseif ~isempty(idCalcReax) && isempty(id_eB),
        solvePathwayID = 1;
    elseif ~isempty(idCalcReax) && ~isempty(id_eB),
        solvePathwayID = 2;
    end
