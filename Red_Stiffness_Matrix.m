% @Author: pradeep
% @Date:   2015-12-04 09:50:53
% @Last Modified by:   pradeep
% @Last Modified time: 2015-12-16 16:12:43
% function [K_red,F_red,U_Known,Disp_Flag] = Red_Stiffness_Matrix(K,F,Left_Edge,Right_Edge,Top_Edge,Bottom_Edge,Node_BR,Node_BL,...
%                                               Node_TL,Node_TR,Ndof,Num_Nodes)
function [K_red,F_red,U_Known,Disp_Flag,U_red_old] = Red_Stiffness_Matrix(K,F,Ndof,Num_Nodes,U_old_eles,N_bnd_mat,bottom_nodes,...
                                                                top_nodes,right_nodes,left_nodes,U_old_inp)

Disp_Flag = zeros(Num_Nodes*Ndof,1);


for j=1:length(bottom_nodes)
  Disp_Flag(bottom_nodes(j)*2) = 1;
end

for j=1:length(left_nodes)
  Disp_Flag( (left_nodes(j)*2)-1) = 1;
end

Ndof_Constraint = sum(Disp_Flag);

U_Overall = zeros(Num_Nodes*Ndof,1);

U_Known = zeros(Ndof_Constraint,1);

count = 1;
for i=1:Num_Nodes*Ndof
  if Disp_Flag(i) == 1
    U_Known(count) = U_Overall(i);
    count = count + 1;
  end
end

F_red = zeros(Num_Nodes*Ndof - Ndof_Constraint,1);
U_red_old = zeros(Num_Nodes*Ndof - Ndof_Constraint,1);
K_red = zeros(Num_Nodes*Ndof - Ndof_Constraint);
K_vec = zeros(Num_Nodes*Ndof,Ndof_Constraint);

count = 1;
count_const = 1;
for i=1:Num_Nodes*Ndof
    if(Disp_Flag(i)==0)
        F_red(count)=F(i);
        U_red_old(count) = U_old_inp(i);
        count = count+1;
      else
        K_vec(:,count_const) = K(:,i);
        count_const = count_const + 1;
    end
end


K_vec_red  = zeros(Num_Nodes*Ndof - Ndof_Constraint,Ndof_Constraint);

for j=1:Ndof_Constraint
  Row_Count = 1;
  for i=1:Num_Nodes*Ndof
    if Disp_Flag(i) == 0
      K_vec_red(Row_Count,j) = K_vec(i,j);
      Row_Count = Row_Count + 1;
    end
  end
end

F_red = F_red - K_vec_red*U_Known;

% Calculating the reduced stiffness matrix
count1=1;
count2=1;
for ii = 1:(Num_Nodes*Ndof)
 if (Disp_Flag(ii)==0)
   for jj = 1:(Num_Nodes*Ndof)
     if (Disp_Flag(jj)==0)
       K_red(count1,count2)=K(ii,jj);
       count2=count2+1;
     end
   end
   count1=count1+1; 
   count2=1;
 end
end