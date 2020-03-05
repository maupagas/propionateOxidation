function [ConcV, prodConc, posConc, ratiospc] = solveRatio(StM, stoV, n_ATP, constConcV, numProtTransloc, Param, LoopData)
  
%     %Identify ratio species
%     idRatiospc  = strcmp(StType, 'r');
%     idRatioSub  = find((stoV < 0) .* idRatiospc); 
%     idRatioProd = find((stoV > 0) .* idRatiospc);

    %Ratio species (Cn/C1)
    id_Cn = LoopData.id_Cn;
    id_C1 = LoopData.id_C1;
    posConc = id_Cn;

    StM(id_Cn) = Param.Max_Conc;
    loopFlag = 2;  %Allows product to go over parameter max
    varConcV = (constConcV == 0);

    %Calculate product concentration
    [ConcV, prodConc] = calcProdConc(Param, StM, stoV, n_ATP, constConcV, varConcV, numProtTransloc, loopFlag);            

%     if posConc == id_C1,
%         sprintf('Product is correct');
%     else
%         sprintf('Product is not the one we are aiming for');
%     end

    %Calculate ratio between species to see if it exceeds the
    %maximum/minimum ratio possible
    ratiospc = prodConc / Param.Max_Conc;

    %Readjust ratios (if product is a higher concentration, decrease
    %substrate concentration to put them within range)
    if prodConc > Param.Max_Conc, 
        prodConc = Param.Max_Conc;
        ConcV(id_Cn) = prodConc / ratiospc;
    end
        ConcV(id_C1) = prodConc;                                  
end