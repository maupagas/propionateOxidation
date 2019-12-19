function [dp, EV] = calc_dp(St, Param)

%Calculate difference of potential to meet the parameter of DG of a proton
ui_Hout  = Param.DG0ft(St.id.H)+ Param.Rth*Param.T * log(St.StM(St.id.H));

dp = (ui_Hout + Param.DG_Prot - Param.DG0ft(St.id.H) - Param.Rth * Param.T * log(St.StM(St.id.H)))/(Param.chrV(St.id.H) * Param.F) + ...
             Param.Rth*Param.T/Param.F * log(St.StM(St.id.H)/St.StM(St.id.Hout));

%Electrical potential
% EV = (dp - (Param.Rth * Param.T / Param.F)*log(St.StM(St.id.H)/St.StM(St.id.Hout))) * St.compPhase;

EV = (dp - (Param.Rth * Param.T / Param.F)*log(St.StM(St.id.H)/St.StM(St.id.Hout))) * St.compPhase * 0;