% @Author: pradeep
% @Date:   2015-12-03 15:21:39
% @Last Modified by:   pradeep
% @Last Modified time: 2015-12-16 16:04:32

function [Kt,K,F] = Stiffness_Matrix(N,dN_dxi,dN_deta,nodes_ele,n_gp,Nodes,Elements,Ndof,Nele,Nodes_ele,...
						 Num_Nodes,C,N_bnd_mat,sigma_active,U_old_Eles,delta_t,kv_bar,epsdot_0,sig_max,...
						 sigma_ij,NI,eps_dot_homo_int,eta_homo_int,eps_homo_int,eps_old_homo_int,...
						 eps_0,eps_1,N_mat,xi_H_eles,F_rhs_old_inp_eles,lambda_s,shift_const,xi_0)

%%%Calculation of Strain displacement matrix B

K = zeros(Num_Nodes*Ndof);	
Kt = zeros(Num_Nodes*Ndof);	
F = zeros(Num_Nodes*Ndof,1);

[gp_vec,w_gp]=Gauss_Points(n_gp);

%%%%Element stiffness matrices
	for e=1:Nele
		sigma_homogenised = sigma_ij(:,e);
		K_ele = zeros(Nodes_ele*Ndof);
		Kt_ele = zeros(Nodes_ele*Ndof);
		F_ele = zeros(Nodes_ele*Ndof,1);
		ID = Elements(e,:);
		x_vec_ele(:,1) = Nodes(ID,1);
		y_vec_ele(:,1) = Nodes(ID,2);
		F_rhs_old_inp_dir = zeros(2,1);

		ID_Dof = Calc_ID(Ndof,ID);
		F_rhs_old_inp_tmp = F_rhs_old_inp_eles(:,e);

		b = 1e-6;
		for j=1:n_gp
			for i=1:n_gp
				if j==1 && i==1
					F_rhs_inp = F_rhs_old_inp_tmp(1:2);
					xi_H_old_tmp = xi_H_eles(1,e);
				elseif j==1 && i==2
					F_rhs_inp = F_rhs_old_inp_tmp(3:4);
					xi_H_old_tmp = xi_H_eles(2,e);
				elseif j==2 && i==1
					F_rhs_inp = F_rhs_old_inp_tmp(7:8);
					xi_H_old_tmp = xi_H_eles(4,e);
				else
					F_rhs_inp = F_rhs_old_inp_tmp(5:6);
					xi_H_old_tmp = xi_H_eles(3,e);
				end
										
				N_vec = N_mat(i*2-1:i*2,:,j);
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

			B_voigt = [ B(1,:); B(2,:); B(3,:) ];
			
			K_ele = K_ele + b*B_voigt'*C*B_voigt*w_gp(i)*w_gp(j)*J;
			F_ele = F_ele - b*B'*sigma_homogenised*w_gp(i)*w_gp(j)*J - xi_H_old_tmp*N_vec'*F_rhs_inp*w_gp(i)*w_gp(j)*J;

			num_2 = kv_bar*eps_dot_homo_int(e)/epsdot_0 + shift_const;
			den_2 = sqrt(num_2^2 + 1);	

			if eps_homo_int(e) < 0
				Kt_ele = Kt_ele + b*(-2*exp(-(eps_homo_int(e)/eps_0)^2)*(eps_homo_int(e)/eps_0)*(eta_homo_int(e)*sig_max + eta_homo_int(e)*sig_max*num_2/den_2)...
								+ exp(-(eps_homo_int(e)/eps_0)^2)*(eta_homo_int(e)*sig_max*kv_bar/epsdot_0)*(1/( (num_2^2 + 1)^(3/2)))*(1/delta_t))*B'*B*w_gp(i)*w_gp(j)*J...
								+ b*B_voigt'*C*B_voigt*w_gp(i)*w_gp(j)*J + xi_H_old_tmp*lambda_s*N_vec'*N_vec*w_gp(i)*w_gp(j)*J;
			else
				Kt_ele = Kt_ele + b*((-2*exp(-(eps_homo_int(e)/eps_0)^2)*(eps_homo_int(e)/eps_0) + 2*(eps_homo_int(e)/eps_1)*(1/eps_1))*(eta_homo_int(e)*sig_max + eta_homo_int(e)*sig_max*num_2/den_2)...
								+ (exp(-(eps_homo_int(e)/eps_0)^2) + (eps_homo_int(e)/eps_1)^2 )*(eta_homo_int(e)*sig_max*kv_bar/epsdot_0)*(1/( (num_2^2 + 1)^(3/2)))*(1/delta_t))*B'*B*w_gp(i)*w_gp(j)*J...
								+ b*B_voigt'*C*B_voigt*w_gp(i)*w_gp(j)*J + xi_H_old_tmp*lambda_s*N_vec'*N_vec*w_gp(i)*w_gp(j)*J;
			end

			end
	
		end

		for Rows = 1:Nodes_ele*Ndof	
			F(ID_Dof(Rows,1)) = F(ID_Dof(Rows,1)) + F_ele(Rows);
			for Cols = 1:Nodes_ele*Ndof
				K(ID_Dof(Rows,1),ID_Dof(Cols,1)) = K(ID_Dof(Rows,1),ID_Dof(Cols,1)) + K_ele(Rows,Cols);
				Kt(ID_Dof(Rows,1),ID_Dof(Cols,1)) = Kt(ID_Dof(Rows,1),ID_Dof(Cols,1)) + Kt_ele(Rows,Cols);
			end
		end

	end