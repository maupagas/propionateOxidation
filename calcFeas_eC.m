%SCRIPT TABULATING e-Carriers and passing them over. 
%SHOULD REPLACE calc_Feasible_eC

function [ratio_eCM, eC_RatioNames, protTranslocFeasComb, Res_eC_Conc, Res_eCFeas_Conc, eC_ConcTable] = calcFeas_eC(Param, St, Reac)

%Ratios of electron carriers allowed
Min_eC_Ratio = Param.Min_Conc / Param.Max_Conc;
Max_eC_Ratio = Param.Max_Conc / Param.Min_Conc;
% Min_eC_Ratio = 1e-5;
% Max_eC_Ratio = 1e5;

num_eCarriers = length(Param.eCarriers);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   PART 1: TABULATE eCarriers CONCENTRATIONS   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%THIS HAS TO BE DONE FOR EVERY PATHWAY 

%defines a range of proton translocated for tabulating NADH and other
%e-Carriers concentrations
NADH_ProtTranslocTabV = -2:1:2;      
eC_ProtTranslocTabV   = NADH_ProtTranslocTabV;

% protTranslocComb = Reac.(char(Pathway(1))).AllProtTranslocComb;

for id_pathway = Param.firstPath2Eval:Param.lastPath2Eval
    reacPath = Reac.Pathway(id_pathway);
    
    %Tabulate values that reoxidate directly to hydrogen
    for i = 1:length(Param.eCarriers)
        eCarrier = Param.eCarriers(i);
        eC_ConcTab = zeros(1, length(Param.eCarriers));
        
        %Identify position of reoxidation reaction and position od the
        %carrier in the concentration matrix
        id_eC = St.id.(char(eCarrier));
        %Reoxidation reaction name
        reoxReac = Param.Reox_reac(strcmp(Param.eCarriers,eCarrier));
        id_eC_Reac = Reac.(char(reacPath)).id.(char(reoxReac));
        
        if id_eC_Reac ~= -1
            
            %Stoichiometric vector
            stoM_ord = Reac.(char(reacPath)).stoM_ord;
            stoV     = stoM_ord(:,id_eC_Reac);
            
            if (strcmp(eCarrier,'NADH') || strcmp(eCarrier,'NADPH') || strcmp(eCarrier,'Fdred') || strcmp(eCarrier,'F420red'))

                for j = 1:length(NADH_ProtTranslocTabV)
                    H_out = NADH_ProtTranslocTabV(j);
                    [eC_ConcTab(j)] = CalcConc_eC(St.StM, stoV, Param, id_eC, H_out);
                end
            end
        else
            eC_ConcTab = -1;
        end
        %Update Structure
        eC_ConcTable.(char(reacPath)).(char(eCarrier)) = eC_ConcTab;
    end
    
    %Store NADH concentrations to calculate reoxidation of other carriers
    if any(strcmp(Param.eCarriers,'NADH')) 
        Conc_NADH = eC_ConcTable.(char(reacPath)).NADH;
        id_NADH = St.id.NADH;
    end
    
    %Tabulate rest of eCarriers that are reoxidized to NADH
    for i = 1:length(Param.eCarriers)
        eCarrier = Param.eCarriers(i);
        %   reacReoxName = Param.Reox_reac(i);
        
        %Identify position of reoxidation reaction and position od the
        %carrier in the concentration matrix
        id_eC = St.id.(char(eCarrier));
        %Reoxidation reaction name
        reoxReac = Param.Reox_reac(strcmp(Param.eCarriers,eCarrier));
        id_eC_Reac = Reac.(char(reacPath)).id.(char(reoxReac));
        
        if id_eC_Reac ~= -1
            
            %Stoichiometric vector
            stoM_ord = Reac.(char(reacPath)).stoM_ord;
            stoV     = stoM_ord(:,id_eC_Reac);
        
            %FADH2 and UQred are reoxidised to NADH
            if strcmp(eCarrier, 'FADH2') || strcmp(eCarrier, 'UQred')
                
                %Preallocate results of concentration of electron carrier
                eC_ConcTab = zeros(length(eC_ProtTranslocTabV), length(NADH_ProtTranslocTabV));
                
                
                for j = 1:length(eC_ProtTranslocTabV)
                    
                    for k = 1:length(NADH_ProtTranslocTabV)
                        
                        % e-Carriers that are not reoxidized directly to H2 need to be recalculated
                        H_out = eC_ProtTranslocTabV(j);
                        St.StM(id_NADH) = Conc_NADH(k);
                        
                        [conc_eC_red] = CalcConc_eC(St.StM, stoV, Param, id_eC, H_out);
                        eC_ConcTab(j,k) = conc_eC_red;
                        
                    end
                end
                
            elseif strcmp(eCarrier, 'NADH') || strcmp(eCarrier, 'Fdred')
                eC_ConcTab = eC_ConcTable.(char(reacPath)).(char(eCarrier));                
            end
        else
            eC_ConcTab = -1;

        end
            %Update Structure
            eC_ConcTable.(char(reacPath)).(char(eCarrier)) = eC_ConcTab;
        
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   PART 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pass up values from tabulated eC (from Part 1) to its correspondent
% translocated protons configuration (each row) for every pathway
for id_pathway = Param.firstPath2Eval:Param.lastPath2Eval
    reacPath = Reac.Pathway(id_pathway);
    %     %Calculates total of proton translocated for the pathway
    %     num_ProtTransloc = length(protTranslocComb.(char(reacPath))(:,1));
    %     %Passes up matrix of proton translocation from structure to variable
    %     protTranslocM     =  protTranslocComb.(char(reacPath));
    
    protTranslocM    = Reac.(char(reacPath)).AllProtTranslocComb;
    num_ProtTransloc = length(protTranslocM(:,1));
    %Preallocates (as 0) the calculation for the matrix ratio for each electron
    %carrier and the number of combinations of proton translocation
    ratioM = zeros(num_ProtTransloc, num_eCarriers);
    id_NADH_reox = Reac.(char(reacPath)).id.NADH_reox;
    
    % Passes up the values of each electron Carrier  (ONE  BY ONE)
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
%         if id_eC_reox ~= -1
%             for rowProtTransloc = 1:num_ProtTransloc
%                 %Finds the value of the proton translocation configuration
%                 %of the particular eCarrier
%                 eC_reax_ProtTransloc = protTranslocM(rowProtTransloc, id_eC_reox);
%                 
%                 %Identifier of the ProtonTranslocation for the tabulated
%                 %Matrix (vectorized is much faster than using find)
%                 id_eC_ProtTranslocTab = (eC_ProtTranslocTabV   == eC_reax_ProtTransloc);   
%                 
%                 %Passes up the value from the tabulated eCarrier to the reduced eCarrier (eC_red)
%                 
%                 %                 %Option A: e-Carrier is FADH2 or Ubiquinone (conserved by oxidising NAD+)
%                 
%                 if (id_FADH2 || id_UQred)
%                     protTransloc_NADH    = protTranslocM(rowProtTransloc, id_NADH_reox);
%                     id_NADH_ProtTransloc = NADH_ProtTranslocTabV == protTransloc_NADH;
%                     eC_redV(rowProtTransloc) = tab_eC_Conc(id_eC_ProtTranslocTab, id_NADH_ProtTransloc);
%                     %Option B: e-Carrier is conserved by converting to H2
%                 else
%                     eC_redV(rowProtTransloc) = tab_eC_Conc(id_eC_ProtTranslocTab);
%                 end
% 
%             end
%         end

        %%% ALTERNATIVE CODE
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
                    id_Prot_eC = eC_ProtTranslocV == eC_ProtTranslocTabV(prot_eC);
                    for prot_eC_NADH = 1:length(eC_ProtTranslocTabV)
                        
                        idProtNADH = protTransloc_NADHV == eC_ProtTranslocTabV(prot_eC_NADH);
                        
                        eC_redV(id_Prot_eC & idProtNADH) = tab_eC_Conc(prot_eC, prot_eC_NADH);
                        
%                         eC_redV(eC_ProtTranslocV == eC_ProtTranslocTabV(prot_eC) & protTransloc_NADHV == eC_ProtTranslocTabV(prot_eC_NADH)) = tab_eC_Conc(prot_eC, prot_eC_NADH);
                    end
                end
                
                %Option B: e-Carrier is conserved by converting to H2
            else
                for prot_eC = 1:length(eC_ProtTranslocTabV)
                    id_Prot_eC = eC_ProtTranslocV == eC_ProtTranslocTabV(prot_eC);
                    eC_redV(id_Prot_eC) = tab_eC_Conc(prot_eC);
%                     eC_redV(eC_ProtTranslocV == eC_ProtTranslocTabV(prot_eC)) = tab_eC_Conc(prot_eC);
                end
            end
            
            
        end
        
        %Pass up the value of the reduced eCarrier and the calculated ratio
        %for only thsoe values that are not 0
        if all(eC_redV)
            Res_eC_Ratio.(char(reacPath)).(char(eCarrier)) = eC_redV / Conc_eC_ox;
            Res_eC_Conc.(char(reacPath)).(char(eCarrier)) = eC_redV;
            ratioM(:,num_eC) = eC_redV / Conc_eC_ox;
        end
        
    end


    %Trim Results
    non_Reac_eC = (ratioM(1,:) == 0);
    ratioM(:,non_Reac_eC) = [];
    
    %Identify feasible results to trim proton translocation combination
    
    % Find eCarriers that are within the feasible range
    idRatios = and(ratioM < Max_eC_Ratio, ratioM > Min_eC_Ratio);
    % Sum all eCarriers that are within the range
    sum_idRatios = sum(idRatios, 2);
    % Identify those combinations of eCarriers where ALL eCarriers are within range
    eC_Feasibility = (sum_idRatios == length(ratioM(1,:)));
    % Collect feasible eCarriers (Only those that have ALL eCarriers feasible)
    ratioFeasM = ratioM(eC_Feasibility == 1, :);
    
    %Pass values to the major structures
    allRatio_eCM.(char(reacPath))         = ratioM;
    %Feasible eCarriers
    ratio_eCM.(char(reacPath))            = ratioFeasM;
    eC_RatioNames.(char(reacPath))        = Param.Reox_reac(non_Reac_eC == 0);    %%%% FIX THIS %%%%%%
    protTranslocFeasComb.(char(reacPath)) = protTranslocM (eC_Feasibility == 1,:);
    eC_Path = fieldnames(Res_eC_Conc.(char(reacPath)));


for n = 1:length(eC_Path)
    eCarrier = eC_Path(n);
    Res_eCFeas_Conc.(char(reacPath)).(char(eCarrier)) = Res_eC_Conc.(char(reacPath)).(char(eCarrier))(eC_Feasibility == 1,:);
end

end