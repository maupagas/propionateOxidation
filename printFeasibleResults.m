function [PrintResults] = printFeasibleResults(Results, Reac, Param, combinationListNames)

%Trial script to print results
% pathwayList = {'F4', 'F5', 'F6'};
% combinationListNames = {'Combination_1', 'Combination_2'};
% netHTransloc = Param.netHTranslocNames;
% combinationListNames(3:end) = [];
printCounter = 0;
PrintResults.FeasComb = [];

% netHTranslocList = {'NetHTransloc__1', 'NetHTransloc_0'};
for k = Param.firstPath2Eval:Param.lastPath2Eval
    reacPath = Reac.Pathway(k);
    pathwayCounter = 0;
    %Restart the matrices to store results for each pathway
    prevProtTranslocRes = [];
    prev_netHTransloc   = [];
    prev_netATP         = [];
    prevDGrV_Res        = [];
    prevConcV_Res       = [];
    
    prevInputsM_Res     = [];
    prev_eC_Ratios_Res  = [];
    prevdp_Res          = [];
    prevDG_Prot_Res     = [];
    prevDGr_FullSto_Res = [];
    
    %Store results for each feasible net proton translocation obtained
    for j = 1:length(combinationListNames)
        
        combinationName = combinationListNames(j);
        
          %Check if there are feasible results for that pathway and combination
        %Pre-write names
        if isempty(Results.(char(combinationName)).(char(reacPath)).CompCombId)
        
        else
        compCombId = Results.(char(combinationName)).(char(reacPath)).CompCombId;
      
%         if ~isempty(compCombId)
            
            %Proton Translocation Combination for Reasible reactions
            protTranslocRes = Results.(char(combinationName)).(char(reacPath)).ProtTranslocComb(compCombId,:);
            
            %Net Proton Translocation
            netHTransloc_Res = Results.(char(combinationName)).(char(reacPath)).netHTransloc(compCombId,:);
            
            %Net ATP Produced
            netATP_Res = Results.(char(combinationName)).(char(reacPath)).netATP_Prod(compCombId,:);
            
            %Gibbs Free Energies
            DGrV_Res = Results.(char(combinationName)).(char(reacPath)).DGrV(compCombId,:);
            
            %Concentrations
            ConcV_Res = Results.(char(combinationName)).(char(reacPath)).OptimConc(:,compCombId)';
            
            %Inputs
            Inputs_Res = Results.(char(combinationName)).(char(reacPath)).Inputs(compCombId,:);
            
            %eC_Ratios
            eC_Ratios_Res =  Results.(char(combinationName)).(char(reacPath)).eC_ratios(compCombId,:);
            
            %Cell potential
            dp_Res = Results.(char(combinationName)).(char(reacPath)).dp(compCombId,:);
            
            %Gibbs of a proton
            DG_Prot_Res = Results.(char(combinationName)).(char(reacPath)).DG_Prot(compCombId,:);
            
            %Gibbs full reaction
            DGr_FullSto_Res = Results.(char(combinationName)).(char(reacPath)).DGr_FullSto(compCombId,:);
            
            %Print all results
            if printCounter == 0 && pathwayCounter == 0
                PrintResults.FeasComb.(char(reacPath)).ProtTranslocComb = protTranslocRes;
                PrintResults.FeasComb.(char(reacPath)).netHTranslocV    = netHTransloc_Res;
                PrintResults.FeasComb.(char(reacPath)).netATPV          = netATP_Res;
                PrintResults.FeasComb.(char(reacPath)).DGrV             = DGrV_Res;
                PrintResults.FeasComb.(char(reacPath)).ConcV            = ConcV_Res;
                
                %Extra information
                PrintResults.FeasComb.(char(reacPath)).Inputs        = Inputs_Res;
                PrintResults.FeasComb.(char(reacPath)).eC_Ratios_Res = eC_Ratios_Res;
                PrintResults.FeasComb.(char(reacPath)).dp            = dp_Res;
                PrintResults.FeasComb.(char(reacPath)).DG_Prot       = DG_Prot_Res;
                PrintResults.FeasComb.(char(reacPath)).DGr_FullSto   = DGr_FullSto_Res;
                
                printCounter = 1;
                pathwayCounter = 1;
                
            else
                PrintResults.FeasComb.(char(reacPath)).ProtTranslocComb = [prevProtTranslocRes; protTranslocRes];
                PrintResults.FeasComb.(char(reacPath)).netHTranslocV    = [prev_netHTransloc;   netHTransloc_Res];
                PrintResults.FeasComb.(char(reacPath)).netATPV          = [prev_netATP;      netATP_Res];
                PrintResults.FeasComb.(char(reacPath)).DGrV             = [prevDGrV_Res;        DGrV_Res];
                PrintResults.FeasComb.(char(reacPath)).ConcV            = [prevConcV_Res;       ConcV_Res];
                
                PrintResults.FeasComb.(char(reacPath)).Inputs        =  [prevInputsM_Res; Inputs_Res] ;
                PrintResults.FeasComb.(char(reacPath)).eC_Ratios_Res =  [prev_eC_Ratios_Res; eC_Ratios_Res];
                PrintResults.FeasComb.(char(reacPath)).dp            =  [prevdp_Res; dp_Res];
                PrintResults.FeasComb.(char(reacPath)).DG_Prot       =  [prevDG_Prot_Res; DG_Prot_Res];
                PrintResults.FeasComb.(char(reacPath)).DGr_FullSto   =  [prevDGr_FullSto_Res; DGr_FullSto_Res];
            end
            
            %Store results to append matrices afterwards
            prevProtTranslocRes = PrintResults.FeasComb.(char(reacPath)).ProtTranslocComb;
            prev_netHTransloc   = PrintResults.FeasComb.(char(reacPath)).netHTranslocV;
            prev_netATP         = PrintResults.FeasComb.(char(reacPath)).netATPV;
            prevDGrV_Res        = PrintResults.FeasComb.(char(reacPath)).DGrV;
            prevConcV_Res       = PrintResults.FeasComb.(char(reacPath)).ConcV;
            
            prevInputsM_Res     = PrintResults.FeasComb.(char(reacPath)).Inputs;
            prev_eC_Ratios_Res  = PrintResults.FeasComb.(char(reacPath)).eC_Ratios_Res;
            prevdp_Res          = PrintResults.FeasComb.(char(reacPath)).dp;
            prevDG_Prot_Res     = PrintResults.FeasComb.(char(reacPath)).DG_Prot;
            prevDGr_FullSto_Res = PrintResults.FeasComb.(char(reacPath)).DGr_FullSto;
            
            
        end
    end
end

end


%Collect total proton translocation
%     if isfield(PrintResults.FeasComb, pathway),
%         %%% THIS DOESN'T HAVE TO BE CALCULATED HERE! CALCULATE NET SOMEWHRE
%         %%% ELSE AND FILTER HERE
%         PrintResults.FeasComb.(char(pathway)).netHTransloc = sum(PrintResults.FeasComb.(char(pathway)).ProtTranslocComb,2) + Param.n_H_ATP;
%     end


