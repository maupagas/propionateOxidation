function [ConcV, prodConc, posConc, ratiospc] = solveRatioFerredoxin(StM, stoV, n_ATP, constConcV, numProtTransloc, Param, St)
  
%     %Identify ratio species
%     idRatiospc  = strcmp(StType, 'r');
%     idRatioSub  = find((stoV < 0) .* idRatiospc); 
%     idRatioProd = find((stoV > 0) .* idRatiospc);

    %Ratio species (Cn/C1)
    id_Fdox  = St.id.Fdox;
    id_Fdred = St.id.Fdred;

    StM(id_Fdox) = Param.Max_Conc;
    loopFlag = 2;  %Allows product to go over parameter max

    %Calculate product concentration
    [ConcV, prodConc, posConc] = calcProdConc(Param, StM, stoV, n_ATP, constConcV, numProtTransloc, loopFlag);            

    %Calculate ratio between species to see if it exceeds the
    %maximum/minimum ratio possible
    ratiospc = prodConc / Param.Max_Conc;

    %Readjust ratios (if product is a higher concentration, decrease
    %substrate concentration to put them within range)
    if prodConc > Param.Max_Conc, 
        prodConc = Param.Max_Conc;
        ConcV(id_Fdox) = prodConc / ratiospc;
    end
        ConcV(id_Fdred) = prodConc;                                  
end