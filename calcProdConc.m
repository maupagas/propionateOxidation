function [concV, prodConc] = calcProdConc(Param, StM, stoV, n_ATP, constConcV, varConcV, numProtTransloc, loopFlag)
           
% Need to redefine this part of the code to load the structures from "Reac" -->  
%Define which variables change and which do not
% constConcV = constConcV == 1;
DG0ft_const  =  Param.DG0ft(constConcV);
DG0ft_var    =  Param.DG0ft(varConcV);
StM_const    = StM(constConcV);
stoV_const   = stoV(constConcV);


%Calculate the Gibbs required 
ObjDG = 0 - numProtTransloc * Param.DG_Prot;
% posConc = find(varConcV == 1);

% If the reaction does not have any product to change (e.g.: Acetate
% Outside of the Cell or e-Carriers concentrations, Calculate as -1 for the
% position and product concentration and DO NOT UPDATE the concentrations
% Vector
% if sum(varConcV) == 0
% %     posConc = -1; 
%     prodConc = -1;
%     concV = StM;
% else
%If there is a product to change, calculate the concentration of the product required for a DG of 0. (ONLY WORKS WITH SINGLE PRODUCTS, NOT WITH MULTIPLE)
%     u_i = DG0ft(constConcV) + Rth * T * log(StM(constConcV)) + chrV(constConcV)  * F .* EV(constConcV);
u_i = DG0ft_const + Param.Rth * Param.T * log(StM_const);

%Sums all the chemical potentials except the species that is going to
%change. Also adds the equivalent Gibbs of the translocated protons
sum_u_i = sum(u_i .* stoV_const) + ObjDG - Param.DG_ATP * n_ATP;


% prodStoich = stoV(varConcV==1) ;        % Check which one of the two is faster
%     prodConc = exp(((-sum_u_i - DG0ft(varConcV) - chrV(varConcV)  * F .* EV(varConcV)) / (Rth * T )));
%Calculates the Product concentration for a DG = 0
prodStoich = sum(stoV .* varConcV); 
prodConc = exp(((-sum_u_i/prodStoich - DG0ft_var) / (Param.Rth * Param.T )));

%Do not allow Product Concentration to go over maximum feasible
%concentration inside the cell.
if prodConc > Param.Max_Conc && loopFlag == 0
    prodConc = Param.Max_Conc;
elseif prodConc > 1 && loopFlag == 1
    %         prodConc = 1;
elseif loopFlag == 2
    %         sprintf('Product concentrations is %.2d', prodConc)
    
end

%Passes up the concentrations to a vector and updates the product
%concentration
%     concV = StM;
StM(varConcV) = prodConc;
concV = StM;




