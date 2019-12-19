function [concV, prodConc, posConc] = calcProdConc(Param, StM, stoV, n_ATP, constConcV, numProtTransloc, loopFlag)

DG0ft = Param.DG0ft;          Rth  = Param.Rth;
T     = Param.T;              chrV = Param.chrV;   
F     = Param.F;               EV = Param.EV;
DG_Prot = Param.DG_Prot;    DG_ATP  = Param.DG_ATP;     maxConc = Param.Max_Conc;

constConcV = constConcV == 1;
varConcV   = constConcV == 0;

%Calculate the Gibbs required 
ObjDG = 0 - numProtTransloc * DG_Prot;
posConc = find(varConcV == 1);

% If the reaction does not have any product to change (e.g.: Acetate
% Outside of the Cell or e-Carriers concentrations, Calculate as -1 for the
% position and product concentration and DO NOT UPDATE the concentrations
% Vector
if sum(varConcV) == 0
    posConc = -1; 
    prodConc = -1;
    concV = StM;
else
    %If there is a product to change, calculate the concentration of the product required for a DG of 0. (ONLY WORKS WITH SINGLE PRODUCTS, NOT WITH MULTIPLE)    
%     u_i = DG0ft(constConcV) + Rth * T * log(StM(constConcV)) + chrV(constConcV)  * F .* EV(constConcV);
    u_i = DG0ft(constConcV) + Rth * T * log(StM(constConcV));

    %Sums all the chemical potentials except the species that is going to
    %change. Also adds the equivalent Gibbs of the translocated protons
    sum_u_i = sum(u_i .* stoV(constConcV)) + ObjDG - DG_ATP * n_ATP;    

    %Calculates the Product concentration for a DG = 0
%     prodConc = exp(((-sum_u_i - DG0ft(varConcV) - chrV(varConcV)  * F .* EV(varConcV)) / (Rth * T )));
    prodConc = exp(((-sum_u_i - DG0ft(varConcV)) / (Rth * T )));

    %Do not allow Product Concentration to go over maximum feasible
    %concentration inside the cell.
    if prodConc > maxConc && loopFlag == 0
        prodConc = maxConc;
    elseif prodConc > 1 && loopFlag == 1
%         prodConc = 1;
    elseif loopFlag == 2
%         sprintf('Product concentrations is %.2d', prodConc)   
       
    end
    
    %Passes up the concentrations to a vector and updates the product
    %concentration
    concV = StM;
    concV(posConc) = prodConc;
end


