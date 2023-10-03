% @Author: pradeep
% @Date:   2015-08-20 16:36:10
% @Last Modified by:   pradeep
% @Last Modified time: 2015-11-25 14:04:49

%%%%This file generates calcium signal based on the concentration of IP3 integrins(S).
%%%%IP3 integrins are obtained by solving the reaction-diffusion equation using the same shape
%%%%functions used for displacement discretisation. strain rate and rate of high affinity integrin
%%%%from the previous time step is used. 

m_s = 1e10;
alpha_S = 1;

kd = 5e-3;
s_0 = 1000e18;
lambda_f = 1.0;
lambda_b = 0.5;

t_fin3 =(index_t)*delta_t;


K_s = zeros(Num_Nodes);	
F_s = zeros(Num_Nodes,1);

[gp_vec,w_gp]=Gauss_Points(n_gp);

	for e=1:Nele
		sigma_homogenised = sigma_ij(:,e);
		K_s_ele = zeros(Nodes_ele);
		Kt_s_ele = zeros(Nodes_ele);
		F_s_ele = zeros(Nodes_ele,1);
		ID = Elements(e,:);
		x_vec_ele(:,1) = Nodes(ID,1);
		y_vec_ele(:,1) = Nodes(ID,2);
		F_rhs_old_inp_dir = zeros(2,1);
		eps_dot_inp = eps_dot_old_inp_C(e);
		ID_Dof_s = Calc_ID(1,ID);
		N_vec_s = zeros(1,4);
		b = 1e-6;
		for j=1:n_gp
			for i=1:n_gp
				if j==1 && i==1
					xi_H_dot_inp = xi_H_dot_new_eles(1,e);
					S_old_inp = S_old_eles(1,e);
				elseif j==1 && i==2
					xi_H_dot_inp = xi_H_dot_new_eles(2,e);
					S_old_inp = S_old_eles(2,e);
				elseif j==2 && i==1
					xi_H_dot_inp = xi_H_dot_new_eles(4,e);
					S_old_inp = S_old_eles(4,e);
				else
					xi_H_dot_inp = xi_H_dot_new_eles(3,e);
					S_old_inp = S_old_eles(3,e);
				end


				N_vec = N_mat(i*2-1:i*2,:,j);
				
				N_vec_s(1) = N_vec(1,1);
				N_vec_s(2) = N_vec(1,3);
				N_vec_s(3) = N_vec(1,5);
				N_vec_s(4) = N_vec(1,7);

				dx_dxi = dN_dxi(i,:,j)*x_vec_ele(:,1);
				dy_dxi = dN_dxi(i,:,j)*y_vec_ele(:,1);
				dx_deta = dN_deta(i,:,j)*x_vec_ele(:,1);
				dy_deta = dN_deta(i,:,j)*y_vec_ele(:,1);

				Jacobian = [dx_dxi dy_dxi;dx_deta dy_deta];
				Jacob_inv = inv(Jacobian);
				J = det(Jacobian);
				
				dN_dX_dy = Jacob_inv*[ dN_dxi(i,:,j); dN_deta(i,:,j)];
				dN_dx = dN_dX_dy(1,:);
				dN_dy = dN_dX_dy(2,:);

				if Ndof == 1
					B(1,:) = dN_dx;
					B(2,:) = dN_dy;
				elseif Ndof ==2
					for k=1:Nodes_ele
						B(1,2*k-1) = dN_dx(k);
						B(2,2*k)   = dN_dy(k);
						B(3,2*k-1) = dN_dy(k);
						B(3,2*k)   = dN_dx(k);
					end
				end
				B_s(1,:) = dN_dx;
				B_s(2,:) = dN_dy;

			B_voigt = [ B(1,:); B(2,:); B(3,:) ];
		
			K_s_ele = K_s_ele + (1/delta_t)*N_vec_s'*N_vec_s*w_gp(i)*w_gp(j)*J + m_s*kT*B_s'*B_s*w_gp(i)*w_gp(j)*J ...
							  + kd*N_vec_s'*N_vec_s*w_gp(i)*w_gp(j)*J + eps_dot_inp*N_vec_s'*N_vec_s*w_gp(i)*w_gp(j)*J;
			F_s_ele = F_s_ele +  (1/delta_t)*S_old_inp*N_vec_s'*w_gp(i)*w_gp(j)*J + (alpha_S/b)*max(0,xi_H_dot_inp)*N_vec_s'*w_gp(i)*w_gp(j)*J;	

			end
	%
		end

		for Rows = 1:Nodes_ele	
			F_s(ID_Dof_s(Rows,1)) = F_s(ID_Dof_s(Rows,1)) + F_s_ele(Rows);
			for Cols = 1:Nodes_ele
				K_s(ID_Dof_s(Rows,1),ID_Dof_s(Cols,1)) = K_s(ID_Dof_s(Rows,1),ID_Dof_s(Cols,1)) + K_s_ele(Rows,Cols);
			end
		end
	end

S_new = sparse(K_s)\sparse(F_s);

optionss = odeset('RelTol',1e-8);
[t_out,c_res] = ode23s(@(t_C,y) myode_c(t_C,y,lambda_f,lambda_b,S_new,s_0),[t t_fin3],c_old,optionss);


S_new_eles = zeros(Nodes_ele,Nele);
for j=1:Nele
	ID = Elements(j,:);
	for i=1:Nodes_ele
		S_new_eles(i,j) = S_new(ID(i));
	end
end

S_old = S_new;
S_old_eles = S_new_eles;

[r,c] = size(c_res);
C_new = c_res(r,:)';
C_eles = zeros(Nodes_ele,Nele);

for j = 1:Nele
	ID = Elements(j,:);
    for i = 1:Nodes_ele
     C_eles(i,j) = C_new(ID(i));
    end
end

Cal = zeros(Nele,1);

for j = 1:Nele
	Cal(j) = mean(C_eles(:,j));
end

c_old = C_new;