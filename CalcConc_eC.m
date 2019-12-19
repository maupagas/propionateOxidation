function [conc_eC_red] = CalcConc_eC(StM, stoV, Param, id_eC, H_out)
%      [eC_ConcTab(j)] = CalcConc_eC(St, Reac, Param, eCarrier, reacPath, H_out);

DG0ft   = Param.DG0ft;      Rth  = Param.Rth;
T       = Param.T;          chrV = Param.chrV;   
F       = Param.F;          EV   = Param.EV;
DG_Prot = Param.DG_Prot;    
% stoM_ord = Reac.(char(reacPath)).stoM_ord;

%Calculate the chemical potentials of all the components except the reduced
%version of the carrier
u_i = DG0ft +  Rth * T * log(StM) + chrV * F  .* EV;
% id_eC = St.id.(char(eC));
% 
% %Reoxidation reaction name
% reoxReac = Param.Reox_reac(strcmp(Param.eCarriers,eC));
% id_eC_Reac = Reac.(char(reacPath)).id.(char(reoxReac));
% 
% %Stoichiometric vector
% stoV     = stoM_ord(:,id_eC_Reac);

%Chemical potential and stoichiometry for all the components besides the e-carrier
u_i_non_eC  = u_i;              u_i_non_eC(id_eC) = [];
stoV_non_eC = stoV;            stoV_non_eC(id_eC) = [];

%Calculate the concentration of the reduced carrier
if stoV(id_eC) == 0
    conc_eC_red = 0;
else
%   ui_eC_red = -(sum(u_i([1:St.id.(char(eC))-1 St.id.(char(eC))+1:end]).*stoV([1:St.id.(char(eC))-1 St.id.(char(eC))+1:end]))- H_out * Param.DG_Prot)/stoV(St.id.(char(eC))); 
    ui_eC_red = -(sum(u_i_non_eC .* stoV_non_eC) - H_out * DG_Prot)/stoV(id_eC);    
    conc_eC_red = exp((ui_eC_red  -  DG0ft(id_eC) - chrV(id_eC)  * F .* EV(id_eC)) / (Rth * T));
end