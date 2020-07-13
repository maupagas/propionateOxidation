%Script to write the values of the DG of the reactions and the activities
%required to find the optimum to run the reactions forward (if possible)
xlsxFile ='PropionateResults.xlsx';
file = strcat(pwd, '\', xlsxFile);
headNameSheet = 'P7a';
targetPath = 'P7a';
lastpathway2Print = find(strcmp(char(targetPath), Reac.Pathway));

%Activates Excel server (I needed to add the file location with its full path)
Excel = actxserver ('Excel.Application');
if ~exist(file,'file')
    ExcelWorkbook = Excel.workbooks.Add;
    ExcelWorkbook.SaveAs(file,1);
    ExcelWorkbook.Close(false);
end
invoke(Excel.Workbooks,'Open',file);

%Column to write DG of each reaction (Will write in IT6 first route from IT, AU, AV... and second in IT7, AU7... etc) 
nCell_Conc = 'AK3';
nCell_ReacNames = 'B2';
nCell_ProtTransloc = 'B3';
nCell_ReacNamesDG = 'T2';
nCell_DG ='T3';

%Print all the rs
Result2Print = Results.Combination_1;
    
    nSheet = headNameSheet;
    reacPath = 'P4a';

    % Clear previous data in the Sheet
    Excel.Worksheets.Item(nSheet).Range('B3:CD30000').ClearContents;  
%     Excel.Worksheets.Item(nSheet).Range('B54:AI11000').ClearContents;  

    %Write calculated Concentrations
    xlswrite1(xlsxFile, Result2Print.(char(reacPath)).OptimConc', nSheet, nCell_Conc);

    %Write the names of the reactions involved in the given pathway
    xlswrite1(xlsxFile, Reac.(char(reacPath)).Names_ord,    nSheet, nCell_ReacNames);
    xlswrite1(xlsxFile, Reac.(char(reacPath)).Names_ord,    nSheet, nCell_ReacNamesDG);

    %Write the results calculated for DG
    xlswrite1(xlsxFile, Result2Print.(char(reacPath)).ProtTranslocComb,  nSheet, nCell_ProtTransloc);
    xlswrite1(xlsxFile, Result2Print.(char(reacPath)).DGrV,  nSheet, nCell_DG);



% for j=1:length(Reac.Names),
%     nCol_DG  = strcat('IT',num2str(j+5));
% %Then run the new xlswrite1 function as many times as needed or in a loop (for example xlswrite1(File,data,location). Then run the following code to close the activex server:
%     xlswrite1(File, Reac.(strcat('DGr',num2str(j))), nSheet, nCol_DG);
% end

invoke(Excel.ActiveWorkbook,'Save');
Excel.Quit
Excel.delete 
clear Excel

clearvars -except Reac Param St Results PrintResults idLoop