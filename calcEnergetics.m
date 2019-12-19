%This is a function to calculate the bioenergetics (DGr) of the provided reactions (from ReacM)

%*************************************************************************************************************
%Some parameters are required:
% a_i : Activities of the different components used in the reactions
%
% Param: Structure containing different parameters:
%      .DG0f : Gibbs energy formation of every component (kJ/mol)
%      .Rth  : Universal gas constant (kJ/(mol K))
%      .T    : Temperature (K)
%      .chrV : Vector of the charge for any of the assumed components to be participating in the reactions
%      .F    : Faraday Constant (kC/mol-e)
%      .EV   : Electric potential betwen the inner and the outer space of the cell
%
% ReacM: Provides the stoichiometry of the considered reactions (Each column corresponds to one reaction)
%*************************************************************************************************************

% function [DGrV] = CalcEnergetics(a_i, Param, ReacM)
function [DGrV] = calcEnergetics(a_i, Param, reacM,  n_ATPV, protTransloc)

DG_Prot = Param.DG_Prot;    DG_ATP = Param.DG_ATP;

% if any(a_i == -1),
%     DGrV = 0;
% else
    %Calculate chemical potentials of the components
    ui = Param.DG0ft + Param.Rth * Param.T * log(a_i);
%     ui = Param.DG0f + Param.Rth * Param.T * log(a_i) + Param.chrV * Param.F .* Param.EV;

    %Calculate the Gibbs of the reactions
%     DGrV_prod  = bsxfun(@times, reacM, ui);
    DGrV_prod = reacM .* ui;
    DGrV_sto      = sum(DGrV_prod);
%     %Identify results with imaginary values (Those which concentrations
%     %were defined as -1)
%     id_imag = find(imag(DGrV1) ~= 0, 1); 
    
    DGrV = DGrV_sto  - DG_Prot * protTransloc - DG_ATP * n_ATPV;
    
%     %Find the first positive DG value (0.01 used to ignore values that are
%     %bigger than 0 but negligible -e.g. 1e-14)
%     posDG_positive = find(DGrV>0.01, 1);
    %Replace values after positive DG as 1000 to discard them from the
    %calculation
%      DGrV
