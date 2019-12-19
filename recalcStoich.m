% This function calculates the reoxidation of electron carriers in the system
% for the provided pathway (FADH2, NADH, Fdred)****** Think about how to
% add electron bifurcation calculation
function [adjustedStoich] = recalcStoich(Reac, St, pathway, Param)

%Pass up values initially
stoM_ord       = Reac.(char(pathway)).stoM_ord;
% reacID_pathway = Reac.(char(pathway)).id;
eCReacID       = Reac.(char(pathway)).id;

%Order should be inverted, so NADH is calculated the last one (Remember
%some eC reoxidize by oxidizing NAD+)
eCarriers   = flipud(Param.eCarriers);
reoxNames   = flipud(Param.Reox_reac);
eC_OxList   = flipud(Param.eC_Cons);

%Evaluate 
for i = 1:length(eCarriers)
    
    eC_red = eCarriers(i);            eC_ox = eC_OxList(i);
    reoxName = reoxNames(i);
    idReac = eCReacID.(char(reoxName));
    
    %Identifies both eCarriers pairs (reduced and oxidised form)
    id_eC_red = St.id.(char(eC_red));  id_eC_ox = St.id.(char(eC_ox));
    
    if idReac == -1,        
        continue
    else
         %eC need to be regenerated in the stoichiometry to ensure the
         %conservation of the moiety
        eC_ox_sto  = -sum(stoM_ord(id_eC_ox, 1:idReac-1), 2);
        eC_red_sto = -sum(stoM_ord(id_eC_red,     1:idReac-1), 2);

        stoM_ord(id_eC_ox, idReac)      = eC_ox_sto;
        stoM_ord(id_eC_red,     idReac) = eC_red_sto;
        
         %eC is reoxidized directly to H2
        if (strcmp(eC_red,'FADH2') || strcmp(eC_red,'UQred')),  
           
            stoM_ord(St.id.NAD,  idReac) = eC_red_sto;
            stoM_ord(St.id.NADH, idReac) = eC_ox_sto;

        else   %eC is reoxidized directly to H2 (stoichiometry of the oxidised form is always paired with hydrogen)
            stoM_ord(St.id.H2,   idReac) = eC_ox_sto;     
        end
    end
end


%Close H2 Mass Balance for all the reoxidations of e-Carriers
for i = 1 : length(stoM_ord(1,:))
    stoM_ord(St.id.H, i) = -(sum(Reac.MassBal.H(1:St.id.H-1)   .* stoM_ord(1:St.id.H-1,i)) + ... 
                                           sum(Reac.MassBal.H(St.id.H+1:end) .* stoM_ord(St.id.H+1:end,i)));
end

adjustedStoich = stoM_ord;
