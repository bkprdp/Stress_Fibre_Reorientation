% @Author: pradeep
% @Date:   2015-12-04 10:37:27
% @Last Modified by:   pradeep
% @Last Modified time: 2015-12-04 11:42:23
function [ID_Dof] = Calc_ID(Ndof,ID)
ID_Dof = zeros(Ndof*length(ID),1);

for i=1:length(ID)
	if Ndof==1
		ID_Dof(i) = ID(i);
	elseif Ndof==2
		ID_Dof(i*Ndof-1,1) = ID(i)*Ndof - 1;
		ID_Dof(i*Ndof,1) = ID(i)*Ndof;
	elseif Ndof ==3
		ID_Dof(i*Ndof-2,1) = ID(i)*Ndof -2;
		ID_Dof(i*Ndof-1,1) = ID(i)*Ndof -1;
		ID_Dof(i*Ndof,1)   = ID(i)*Ndof;
	end
end