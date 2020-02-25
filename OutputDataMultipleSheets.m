%Script to write the values of the DG of the reactions and the activities
%required to find the optimum to run the reactions forward (if possible)

%It only can print One combination at a time, otherwise, too many tabs
%would be written

%%%%%%%%%% EDITING OF FILE STARTS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Excel file to print and its path
xlsxFile ='ButyrateResults.xlsx';
file = strcat(pwd, '\', xlsxFile);

printType = 'all';   %choose between 'all' or 'final'
printCombName = 'Combination_18';  %Select Combination List Name

%%%%%%%%%% EDITING OF FILE ENDS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Write the last two missing headers for the excel sheet
nameProtonTransloc = 'Net Proton Translocation';
nameATP_Prod = 'Net ATP Produced';

%Activates Excel server to print in multiple files
Excel = actxserver ('Excel.Application');
if ~exist(file,'file')
    ExcelWorkbook = Excel.workbooks.Add;
    ExcelWorkbook.SaveAs(file,1);
    ExcelWorkbook.Close(false);
end
invoke(Excel.Workbooks,'Open',file);

%Column to write DG of each reaction (Will write in IT6 first route from IT, AU, AV... and second in IT7, AU7... etc) 

%Print results
for i = Param.firstPath2Eval : Param.lastPath2Eval
    
    %Writes in the sheet of the pathway to be printed
    nameSheet = Reac.Pathway(i);
    
    %Print Reaction Names
    nCell_ReacNames    = 'B2';
    ncol_NetHTransloc  = xlscol(nCell_ReacNames) + length(Reac.(char(Reac.Pathway(i))).Names_ord);
    nCell_NetHTransloc = strcat(xlscol(ncol_NetHTransloc), nCell_ReacNames(2));
    nCell_NetATPName   = strcat(xlscol(ncol_NetHTransloc + 1), nCell_ReacNames(2));
    nCell_ReacNamesDG  = strcat(xlscol(ncol_NetHTransloc + 3), nCell_ReacNames(2));

    numColConcNames = xlscol(nCell_ReacNamesDG) + length(Reac.(char(Reac.Pathway(i))).stoM_ord(1,:)) + 1;
    nCell_ConcNames = strcat(xlscol(numColConcNames), nCell_ReacNames(2)); 
         
    numColInputsNames     = numColConcNames + length(St.StM) + 1;
    nCell_InputsNames     = strcat(xlscol(numColInputsNames), nCell_ReacNames(2));
    
    numColeCarrierNames   = numColInputsNames + length(St.combinationListNames) + 1;
    nCell_eCarrierNames   = strcat(xlscol(numColeCarrierNames), nCell_ReacNames(2));
    
    numDGr_FullStoName    = numColeCarrierNames + length(Results.Combination_1.(char(Reac.Pathway(i))).eC_ReoxNames) + 1;
    nCell_DGr_FullStoName = strcat(xlscol(numDGr_FullStoName), nCell_ReacNames(2));
  
    numdp_label           = numDGr_FullStoName + 1;
    nCell_dp_Label        = strcat(xlscol(numdp_label), nCell_ReacNames(2));
    
    num_DG_Prot           = sum(numdp_label) + 1;
    nCell_DG_Prot_Label   = strcat(xlscol(num_DG_Prot),  nCell_ReacNames(2));

    %Print values Names
    nCell_ProtTransloc    = strcat(nCell_ReacNames(1),          nCell_ReacNames(end)       + 1);
    nCell_NetProtTransloc = strcat(nCell_NetHTransloc(1),       nCell_NetHTransloc(end)    + 1);
    nCell_NetATP          = strcat(nCell_NetATPName(1),         nCell_NetATPName(end)      + 1);
    nCell_DG              = strcat(nCell_ReacNamesDG(1),        nCell_ReacNamesDG(end)     + 1);
    nCell_Conc            = strcat(xlscol(numColConcNames),     nCell_ReacNames(end)       + 1);
    nCell_Inputs          = strcat(xlscol(numColInputsNames),   nCell_InputsNames(end)     + 1);
    nCell_eCarriers       = strcat(xlscol(numColeCarrierNames), nCell_eCarrierNames(end)   + 1);
    nCell_DGr_FullSto     = strcat(xlscol(numDGr_FullStoName),  nCell_DGr_FullStoName(end) + 1);
    nCell_dp              = strcat(xlscol(numdp_label),         nCell_dp_Label(end)        + 1);
    nCell_DG_Prot         = strcat(xlscol(num_DG_Prot),         nCell_DG_Prot_Label(end)   + 1);
   
    %Print all the rs
    nSheet  = char(nameSheet);
    pathway = Reac.Pathway(i); 

    % Clear previous data in the Sheet
%     Excel.Worksheets.Item(nSheet).Range('B2:CD30000').ClearContents;  
%     Excel.Worksheets.Item(nSheet).Range('B54:AI11000').ClearContents; 

        switch printType
            case 'all'
                print_ConcResults      = Results.(char(printCombName)).(char(pathway)).OptimConc';
                print_ProtTranslocComb = Results.(char(printCombName)).(char(pathway)).ProtTranslocComb;
                print_DGrV             = Results.(char(printCombName)).(char(pathway)).DGrV;
                
                print_NetProtTransloc  = Results.(char(printCombName)).(char(pathway)).netHTransloc;  
                print_NetATP           = Results.(char(printCombName)).(char(pathway)).netATP_Prod;
                
                print_Inputs           = Results.(char(printCombName)).(char(pathway)).Inputs;
                print_eCarrierRatios   = Results.(char(printCombName)).(char(pathway)).eC_ratios;     %%%%%% FIX THIS ONE
                print_DGr_FullSto      = Results.(char(printCombName)).(char(pathway)).DGr_FullSto;
                print_dp               = Results.(char(printCombName)).(char(pathway)).dp;
                print_DG_Prot          = Results.(char(printCombName)).(char(pathway)).DG_Prot;
                print_eC_Names         = Results.(char(printCombName)).(char(pathway)).eC_ReoxNames;
            
            case 'final'
                               
                if isfield(PrintResults.FeasComb, pathway)
                    print_ConcResults      = PrintResults.FeasComb.(char(pathway)).ConcV;
                    print_ProtTranslocComb = PrintResults.FeasComb.(char(pathway)).ProtTranslocComb;
                    print_DGrV             = PrintResults.FeasComb.(char(pathway)).DGrV;  
                    print_NetProtTransloc  = PrintResults.FeasComb.(char(pathway)).netHTranslocV;  
                    print_NetATP           = PrintResults.FeasComb.(char(pathway)).netATPV;  

                    print_Inputs           = PrintResults.FeasComb.(char(pathway)).Inputs;
                    print_eCarrierRatios   = PrintResults.FeasComb.(char(pathway)).eC_Ratios_Res;     %%%%%% FIX THIS ONE
                    print_DGr_FullSto      = PrintResults.FeasComb.(char(pathway)).DGr_FullSto;
                    print_dp               = PrintResults.FeasComb.(char(pathway)).dp;
                    print_DG_Prot          = PrintResults.FeasComb.(char(pathway)).DG_Prot;
                    print_eC_Names         = Results.(char(printCombName)).(char(pathway)).eC_ReoxNames;

                else
                    print_ConcResults      = [];
                    print_ProtTranslocComb = [];
                    print_DGrV             = [];   
                    print_NetProtTransloc  = [];
                    print_NetATP           = [];
                    print_Inputs           = [];
                    print_eCarrierRatios   = [];
                    print_DGr_FullSto      = [];
                    print_dp               = [];
                    print_DG_Prot          = [];
                    print_eC_Names = [];
                end
        end
        
        %Clean Sheet
        emptyStrings = repmat({''},2000,300);
        xlswrite1(xlsxFile, emptyStrings, nSheet, nCell_ReacNames);

        if isempty(print_ConcResults) == 0
        
            %Write calculated Concentrations
            xlswrite1(xlsxFile, print_ConcResults, nSheet, nCell_Conc);
            xlswrite1(xlsxFile, St.StNames', nSheet, nCell_ConcNames);

            %Write the names of the reactions involved in the given pathway
            xlswrite1(xlsxFile, Reac.(char(pathway)).Names_ord,    nSheet, nCell_ReacNames);
            xlswrite1(xlsxFile, cellstr(nameProtonTransloc),       nSheet, nCell_NetHTransloc);        
            xlswrite1(xlsxFile, Reac.(char(pathway)).Names_ord,    nSheet, nCell_ReacNamesDG);
            xlswrite1(xlsxFile, cellstr(nameATP_Prod),             nSheet, nCell_NetATPName);   

            %Write the results calculated for DG
            xlswrite1(xlsxFile, print_ProtTranslocComb,  nSheet, nCell_ProtTransloc);
            xlswrite1(xlsxFile, print_NetProtTransloc,   nSheet, nCell_NetProtTransloc);
            xlswrite1(xlsxFile, print_NetATP,   nSheet, nCell_NetATP);
            xlswrite1(xlsxFile, print_DGrV,  nSheet, nCell_DG);
            
            %Print Extra Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            xlswrite1(xlsxFile, St.combinationListNames,  nSheet, nCell_InputsNames);           
            xlswrite1(xlsxFile, print_Inputs,   nSheet, nCell_Inputs);
                        
            xlswrite1(xlsxFile, print_eC_Names',        nSheet, nCell_eCarrierNames);
            xlswrite1(xlsxFile, print_eCarrierRatios,  nSheet, nCell_eCarriers);
            
            xlswrite1(xlsxFile, cellstr('DGr'),     nSheet, nCell_DGr_FullStoName);         
            xlswrite1(xlsxFile, print_DGr_FullSto,  nSheet, nCell_DGr_FullSto);
                        
            xlswrite1(xlsxFile, cellstr('dp'),  nSheet, nCell_dp_Label);
            xlswrite1(xlsxFile, print_dp,  nSheet, nCell_dp);
                       
            xlswrite1(xlsxFile, cellstr('DG_Prot'),  nSheet, nCell_DG_Prot_Label);
            xlswrite1(xlsxFile, print_DG_Prot,  nSheet, nCell_DG_Prot);
        
        end

end


%Saves excel datasheet to be printed and closes the XServer
invoke(Excel.ActiveWorkbook,'Save');
Excel.Quit
Excel.delete 
clear Excel

%Clears all the variables except the interested ones
clearvars -except Reac Param idLoop St Results PrintResults