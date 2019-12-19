load calc_eCarriers.mat

%For each 
for num_eC = 1:num_eCarriers
    eCarrier = Param.eCarriers(num_eC);
    
    %Identify if the electron carrier is FADH2 or Ubiquinone (for later
    %use)
    id_FADH2 = strcmp(eCarrier,'FADH2');
    id_UQred = strcmp(eCarrier,'UQred');
    %Recalls the tabulated values for the electron carrier (Done in
    %part 1)
    tab_eC_Conc = eC_ConcTable.(char(reacPath)).(char(eCarrier));
    
    %Create a flag that tags if there is an electron carrier
    if tab_eC_Conc ~= -1
        isTab_eC_Conc = 1;
    else
        isTab_eC_Conc = 0;
    end
    
    %Identifies the oxidised version of the carrier to store its
    %concentration (to calculate the ratio reduced/oxidised later)
    id_eCox = St.id.(char(Param.eC_Cons(num_eC)));
    Conc_eC_ox = St.StM(id_eCox);
    
    %Preallocate results sizes for the electron carrier (is avector)
    eC_redV = zeros(num_ProtTransloc, 1);
    
    %Identifies reoxidation/conservation moiety reaction
    reacReoxName = Param.Reox_reac(num_eC);
    id_eC_reox   = Reac.(char(reacPath)).id.(char(reacReoxName));
    
    %Passes up for EACH proton translocation configuration the value of
    %the tabulated e-Carrier (THIS IS INEFFICIENT, BETTER TO VECTORIZE
    %THE CODE, e.g.: A(eC_reac_protTranslocV == -1) = AAvalue;
    if id_eC_reox ~= -1
        
         eC_ProtTranslocV = protTranslocM(:, id_eC_reox);
         
         %Initialize vector for concentrations of electron carriers 
         eC_redV = 0 .* eC_ProtTranslocV;
         
         %If carrier is FADH2 or UQred, they depend on the concentration of
         %NADH.....
         if (id_FADH2 || id_UQred)
             protTransloc_NADHV    = protTranslocM(:, id_NADH_reox);
             %For each value of the protons translocated and the
             %"translocated state of NADH", use the required concentrations
             for prot_eC = 1:length(eC_ProtTranslocTabV)
                 
                 for prot_eC_NADH = 1:length(eC_ProtTranslocTabV)

                     eC_redV(eC_ProtTranslocV == eC_ProtTranslocTabV(prot_eC) & protTransloc_NADHV == eC_ProtTranslocTabV(prot_eC_NADH)) = tab_eC_Conc(prot_eC, prot_eC_NADH);
                 end
             end

             %Option B: e-Carrier is conserved by converting to H2
         else
             for prot_eC = 1:length(eC_ProtTranslocTabV)
                 eC_redV(eC_ProtTranslocV == eC_ProtTranslocTabV(prot_eC)) = tab_eC_Conc(prot_eC);
             end
         end
         

    end
end