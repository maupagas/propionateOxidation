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
    pathway = Reac.Pathway(k);
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
    for j = 1:length(Param.netHTranslocNames)

        netHTransloc =  Param.netHTranslocNames(j);

        if isfield(Results, netHTransloc)

            %Store results for each combination of concentrations used
            for i = 1:length(combinationListNames)
                
                combinationName = combinationListNames(i);
                
                if ~isempty(Results.(char(combinationName)).(char(pathway)).ProtTranslocComb)

                    %Check if there are feasible results for that pathway and combination
                    %Pre-write names
                    
                    if isempty(Results.(char(combinationName)).(char(pathway)).CompCombId) == 0

                        %Proton Translocation
                        protTranslocComb = Results.(char(combinationName)).(char(pathway)).ProtTranslocComb;


                        %Check if the combination is in the current proton
                        %translocation combination (CHECK POSSIBILITY TO
                        %DECLARE PREVIOUSLY VARIABLES TO PRINT AND MAKE
                        %CODE MUCH SHORTER)
                        if isfield(Results.(char(netHTransloc)), pathway)

                            if isfield(Results.(char(netHTransloc)).(char(pathway)), combinationName)
                                %Proton Translocation
                                feasHTransloc = Results.(char(netHTransloc)).(char(pathway)).(char(combinationName));
                                protTranslocRes = protTranslocComb(feasHTransloc, :);
                                
                                %Net Proton Translocation
                                netHTranslocation = Results.(char(combinationName)).(char(pathway)).netHTransloc;
                                netHTransloc_Res  = netHTranslocation(feasHTransloc, :);
                                
                                %Net ATP Produced
                                netATP = Results.(char(combinationName)).(char(pathway)).netATP_Prod;
                                netATP_Res  = netATP(feasHTransloc, :);

                                %Gibbs Free Energies
                                DGrV = Results.(char(combinationName)).(char(pathway)).DGrV;
                                DGrV_Res = DGrV(feasHTransloc, :);

                                %Concentrations
                                ConcV = Results.(char(combinationName)).(char(pathway)).OptimConc';
                                ConcV_Res = ConcV(feasHTransloc, :);
    
                                %Inputs             
                                Inputs = Results.(char(combinationName)).(char(pathway)).Inputs;
                                Inputs_Res = Inputs(feasHTransloc, :);
                                
                                %eC_Ratios
                                eC_Ratios =  Results.(char(combinationName)).(char(pathway)).eC_ratios;
                                eC_Ratios_Res = eC_Ratios(feasHTransloc, :);
                                
                                %Cell potential
                                dp = Results.(char(combinationName)).(char(pathway)).dp;
                                dp_Res = dp(feasHTransloc, :);
                                
                                %Gibbs of a proton
                                DG_Prot = Results.(char(combinationName)).(char(pathway)).DG_Prot;
                                DG_Prot_Res = DG_Prot(feasHTransloc, :);
                                
                                %Gibbs full reaction
                                DGr_FullSto = Results.(char(combinationName)).(char(pathway)).DGr_FullSto;
                                DGr_FullSto_Res = DGr_FullSto(feasHTransloc, :);

                                %Print all results
                                if printCounter == 0 && pathwayCounter == 0
                                    PrintResults.FeasComb.(char(pathway)).ProtTranslocComb = protTranslocRes;
                                    PrintResults.FeasComb.(char(pathway)).netHTranslocV    = netHTransloc_Res;
                                    PrintResults.FeasComb.(char(pathway)).netATPV          = netATP_Res;
                                    PrintResults.FeasComb.(char(pathway)).DGrV             = DGrV_Res;
                                    PrintResults.FeasComb.(char(pathway)).ConcV            = ConcV_Res;
                                    
                                    %Extra information
                                    PrintResults.FeasComb.(char(pathway)).Inputs        = Inputs_Res;
                                    PrintResults.FeasComb.(char(pathway)).eC_Ratios_Res = eC_Ratios_Res;
                                    PrintResults.FeasComb.(char(pathway)).dp            = dp_Res;
                                    PrintResults.FeasComb.(char(pathway)).DG_Prot       = DG_Prot_Res;
                                    PrintResults.FeasComb.(char(pathway)).DGr_FullSto   = DGr_FullSto_Res;                         
                                 
                                    printCounter = 1;
                                    pathwayCounter = 1;

                                else
                                    PrintResults.FeasComb.(char(pathway)).ProtTranslocComb = [prevProtTranslocRes; protTranslocRes]; 
                                    PrintResults.FeasComb.(char(pathway)).netHTranslocV    = [prev_netHTransloc;   netHTransloc_Res];
                                    PrintResults.FeasComb.(char(pathway)).netATPV          = [prev_netATP;      netATP_Res];
                                    PrintResults.FeasComb.(char(pathway)).DGrV             = [prevDGrV_Res;        DGrV_Res]; 
                                    PrintResults.FeasComb.(char(pathway)).ConcV            = [prevConcV_Res;       ConcV_Res];   
                                    
                                    PrintResults.FeasComb.(char(pathway)).Inputs        =  [prevInputsM_Res; Inputs_Res] ;
                                    PrintResults.FeasComb.(char(pathway)).eC_Ratios_Res =  [prev_eC_Ratios_Res; eC_Ratios_Res];                                    
                                    PrintResults.FeasComb.(char(pathway)).dp            =  [prevdp_Res; dp_Res];
                                    PrintResults.FeasComb.(char(pathway)).DG_Prot       =  [prevDG_Prot_Res; DG_Prot_Res];
                                    PrintResults.FeasComb.(char(pathway)).DGr_FullSto   =  [prevDGr_FullSto_Res; DGr_FullSto_Res];
                                end

                                    %Store results to append matrices afterwards
                                    prevProtTranslocRes = PrintResults.FeasComb.(char(pathway)).ProtTranslocComb;
                                    prev_netHTransloc   = PrintResults.FeasComb.(char(pathway)).netHTranslocV;  
                                    prev_netATP         = PrintResults.FeasComb.(char(pathway)).netATPV; 
                                    prevDGrV_Res        = PrintResults.FeasComb.(char(pathway)).DGrV;
                                    prevConcV_Res       = PrintResults.FeasComb.(char(pathway)).ConcV;
                                    
                                    prevInputsM_Res     = PrintResults.FeasComb.(char(pathway)).Inputs;
                                    prev_eC_Ratios_Res  = PrintResults.FeasComb.(char(pathway)).eC_Ratios_Res;
                                    prevdp_Res          = PrintResults.FeasComb.(char(pathway)).dp;
                                    prevDG_Prot_Res     = PrintResults.FeasComb.(char(pathway)).DG_Prot;
                                    prevDGr_FullSto_Res = PrintResults.FeasComb.(char(pathway)).DGr_FullSto;
                                    
                                    
                            end
                        end
                    end
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
        
end
