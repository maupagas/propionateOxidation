function [DG_diss] = MinDG_diss(a_i, St, Param, Reac, ReacV)

%Updates the activity values (to ensure that only the desired variables are 
%changed during the use of fmincon
% St.StM(1:Reac.id.H2-1) = a_i(1:Reac.id.H2-1);
a_i = a_i .* (Param.Conc2Change~= 0)+ St.StM .* (Param.Conc2Change==0);
% a_i = St.StM;
%find positions for different activities
% Pos_Var_ai = find(Param.Conc2Change==1);
% for i = 1:length(Pos_Var_ai), 
%     a_i(Pos_Var_ai(i)) = Var_ai(i); 
% end

ReacM = ReacV;

% Mode  = strcat('Mode_', Reac.Modes(i));
% ReacM =  Reac.stoM.(char(Mode))(:,1);

%Calculate the Gibbs of the reactions
[DGrV] = CalcEnergetics(a_i, Param, ReacM);
 
% %Minimize DG_diss
DG_diss = sum(abs(DGrV));

end