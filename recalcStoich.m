% This function calculates the reoxidation of electron carriers in the system
% for the provided pathway (FADH2, NADH, Fdred)****** Think about how to
% add electron bifurcation calculation
function [adjustedStoich] = recalcStoich(Reac, St, pathway, Param)

%Pass up values initially
stoM_ord   = Reac.(char(pathway)).stoM_ord;
eCReacID   = Reac.(char(pathway)).id;
Labels_ord = Reac.(char(pathway)).Labels_ord;
%% Recalculate total products being transported
%Recalculate Final Product Reaction
idTr = find(strcmp(Reac.(char(pathway)).Labels_ord, 'Fp'));
TrR  = stoM_ord(:,idTr);

%Auxiliary multiplicative vector
auxV = 0*TrR;
auxV(1:St.id.H2-1) = 1;

%Substrate identifier
idSub = find(double(TrR < 0) .* auxV);

% Find the final stoichiometry that needs to be transported outside the
% cell 
%(IN THE EXCEL THE TRANSPORT REACTION NEEDS TO BE INITIALLY SET AT 1:1 STOICHIOMETRY 
% TO BE ABLE TO APPLY CORRECTLY THE MULTIPLIER!
stoTrVal = sum(stoM_ord(idSub,1:idTr-1));
stoM_ord(:,idTr) = TrR * stoTrVal;


%% Recalculate e-Carriers Stoichiometry

%Order should be inverted, so NADH is calculated the last one (Remember
%some eC reoxidize by oxidizing NAD+)
eCarriers   = flipud(Param.eCarriers);
reoxNames   = flipud(Param.Reox_reac);
eC_OxList   = flipud(Param.eC_Cons);

% Recalculate electron carriers
for i = 1:length(eCarriers)
    
    eC_red = eCarriers(i);            eC_ox = eC_OxList(i);
    reoxName = reoxNames(i);
    idReac = eCReacID.(char(reoxName));
    
    %Identifies both eCarriers pairs (reduced and oxidised form)
    id_eC_red = St.id.(char(eC_red));  id_eC_ox = St.id.(char(eC_ox));
    
    if idReac == -1        
        continue
    else
         %eC need to be regenerated in the stoichiometry to ensure the
         %conservation of the moiety
        eC_ox_sto  = -sum(stoM_ord(id_eC_ox, 1:idReac-1), 2);
        eC_red_sto = -sum(stoM_ord(id_eC_red,     1:idReac-1), 2);

        stoM_ord(id_eC_ox,  idReac) = eC_ox_sto;
        stoM_ord(id_eC_red, idReac) = eC_red_sto;
        
        %If the reaction has a swap with NADH (eCS), adjust stoichiometry to NADH
        if strcmp(Labels_ord(idReac),'eCS')          
%         if strcmp(eC_red,'FADH2')          
            stoM_ord(St.id.NAD,  idReac) = eC_red_sto;
            stoM_ord(St.id.NADH, idReac) = eC_ox_sto;

        else   %If Label is eC, then the reozidation is directly to H2 (stoichiometry of the oxidised form is always paired with hydrogen)
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
