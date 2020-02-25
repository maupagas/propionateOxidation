function [ID] = idLoops(St, Reac) 

numPathways = length(Reac.Pathway);

for k = 1:numPathways

    reacPath = Reac.Pathway(k);
    stoM = Reac.(char(reacPath)).stoM_ord;
    
    n_ATPV       = Reac.(char(reacPath)).n_ATP_ord;
    protTransloc = Reac.(char(reacPath)).AllProtTranslocComb;
    constConcM   = Reac.(char(reacPath)).constConcM;
    varConcM     = Reac.(char(reacPath)).varConcM;
    posConcV     = Reac.(char(reacPath)).posConcV;
    
    %Substrates and Products Identifiers
    idSubs  = stoM(1:St.id.H2,:) < 0;
    idProds = stoM(1:St.id.H2,:) > 0;

    numSubs  = sum(idSubs , 1);
    numProds = sum(idProds, 1); 
    
    %Identify Loops and ratio calculations
    idRatioRc = strcmp(Reac.(char(reacPath)).Labels_ord, 'Cr');
    idLoops   = strcmp(Reac.(char(reacPath)).Labels_ord, 'Lc');
    
    idLoops = idRatioRc + idLoops;

    
    %Possibility to discount here the substrates and products defined as a
    %ratio

    %Identify number of loops in the stoichiometryid
    numLoops = sum(idLoops); 
    ID.(char(reacPath)).numLoops   = numLoops;
      
    %% Loops identifications
    startPointLoop = find(idLoops);
    LoopsNames = cell(numLoops, 1);

    %Identify start and end of each loop (for pathways that contain a loop)
    if numLoops > 0
        for i = 1:numLoops

            LoopsNames{i} = sprintf('Loop_%d',i);
            startLoop = startPointLoop(i);
            %Identify Start Of Loop
            ID.(char(reacPath)).(char(LoopsNames{i})).Start = startLoop;

            %Find both substrates for the start of the loop
            if idRatioRc(startLoop) == 1
                endLoop = startLoop + 1;
                ID.(char(reacPath)).(char(LoopsNames{i})).TypeLoop = 'ratio';

            else
                ID.(char(reacPath)).(char(LoopsNames{i})).TypeLoop = 'loop';
                loopSubs  = idSubs(:,startLoop);
                loopProds = idProds(:,startLoop);
                %Find the column where the loop ends from its start 
                %(the reaction before Cout a substrate at some point)
    %             [~, endcolLoop] = find(idProds(loopSubs, startLoop:end));
    %             endLoop = endcolLoop + startLoop - 1;

%                 endLoop = find(stoM(loopProds, :) == -1) - 1;
                [~, loopProdCons] = find(stoM(loopProds,:) < -eps);
                endLoop = max(loopProdCons) -1;
%                 endLoop = endLoop(endLoop > startLoop); 
            end
            
            ID.(char(reacPath)).(char(LoopsNames{i})).End   = endLoop;
            stoLoopM = stoM(:, startLoop : endLoop);
            
            %Find in the loop where is Cin, Cout, Cn and C1
            idSubsLoop  = idSubs (:, startLoop);
            idProdsLoop = idProds(:, startLoop);
            %Substrates
            id_Cin = (sum(stoLoopM(1:St.id.H2,:), 2) .* idSubsLoop) <0;
            id_Cn = (idSubsLoop .* (id_Cin == 0)) == 1;
            %Products
            id_Cout = sum(stoLoopM(1:St.id.H2,:), 2) .* idProdsLoop >0;
            id_C1 = (idProdsLoop .* (id_Cout == 0)) == 1;
            
            %Build loop structure matrices and values
            %Matrices
            ID.(char(reacPath)).(char(LoopsNames{i})).stoLoopM          = stoLoopM;
            ID.(char(reacPath)).(char(LoopsNames{i})).n_ATPLoopV        = n_ATPV      (:, startLoop:endLoop);
            ID.(char(reacPath)).(char(LoopsNames{i})).protTranslocLoopV = protTransloc(:, startLoop:endLoop);
            ID.(char(reacPath)).(char(LoopsNames{i})).constConcLoopM    = constConcM  (:, startLoop:endLoop); 
            ID.(char(reacPath)).(char(LoopsNames{i})).varConcLoopM      = varConcM(:, startLoop:endLoop);
            ID.(char(reacPath)).(char(LoopsNames{i})).posConcLoopV     = posConcV(startLoop:endLoop);
            
            %Values
            ID.(char(reacPath)).(char(LoopsNames{i})).id_Cin  = find(id_Cin);
            ID.(char(reacPath)).(char(LoopsNames{i})).id_Cn   = find(id_Cn);
            ID.(char(reacPath)).(char(LoopsNames{i})).id_C1   = find(id_C1);
            ID.(char(reacPath)).(char(LoopsNames{i})).id_Cout = find(id_Cout);                 
        end
    end
end
