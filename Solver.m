%%%%Square cell contracts due to the growth of stress fibres, solved using staggered solver 

clear all
close all
clc
warning on

Ndof = 2;
n_gp = 2;
nu = 0.3;
kv_bar = 6;
sig_max = 20.0e3;
E = 0.08e3;
epsdot_0 = 2.8e-4;
theta = 100;

[Nodes,Elements,Nele,Nodes_ele,Num_Nodes,bottom_nodes,top_nodes,right_nodes,left_nodes] = mesh_file	;

phi = [-pi/2:pi/10:pi/2];
eta_old = zeros(Nele,1);
sigma_active = zeros(Nele,length(phi));
Eps_avg_old = zeros(3,Nele);
U_old_Eles = zeros(Nodes_ele*Ndof,Nele);
sigma_ij = zeros(3,Nele);
eta_phi = zeros(Nele,length(phi));
eta_phi_old = zeros(Nele,length(phi));
eps_dot_homo = zeros(Nele,length(phi));
eps_homo= zeros(Nele,length(phi));
eps_old_homo= zeros(Nele,length(phi));
sigma_active= zeros(Nele,length(phi));
NI = zeros(3,Nele);
U = zeros(Num_Nodes*Ndof,1);
U_old_inp = zeros(Num_Nodes*Ndof,1);
U_old = zeros(Num_Nodes*Ndof,1);
eta_int = zeros(Nele,1);
eps_dot_homo_int= zeros(Nele,1);
eps_homo_int= zeros(Nele,1);
eps_old_homo_int= zeros(Nele,1);
eta_homo_int = zeros(Nele,1);
stretch_delta_old = zeros(Ndof*Num_Nodes,1);

t = 0;
t_fin = 1001;
delta_t = 1;
index_t = 1;
eps_0 = 1e-1;
eps_1 = 1.4*eps_0;
shift_const = 0;

stretch_delta_max = 130e-9;
lambda_s = 0.015e-3;

boltzman_const = 1.38*10^-23;
Temp = 310;                            %Absolute temperature

kT = boltzman_const*Temp;
delta_mu = 5*kT;   
xi_0 = 3000e12;                        %Total concentration of integrins at equilibrium
xi_H_0 = xi_0/(1 + exp(delta_mu/kT));  %Initial concentration of high affinity binders
xi_l_0 = xi_0/(1 + exp(-delta_mu/kT)); %Initial concentration of low affinity binders

xi_H_new = xi_H_0*ones(Num_Nodes,1);
xi_H_old = xi_H_0*ones(Num_Nodes,1);
xi_H_eles = xi_H_0*ones(Nodes_ele*Ndof,Nele);
F_rhs_old_inp_dir = zeros(Ndof*Nodes_ele,Nele);
F_rhs_old_inp_eles = zeros(Ndof*Nodes_ele,Nele);

xi_H_dot_new_eles = zeros(Nodes_ele,Nele);
s_0 = 1000e18;
S_old_eles = ones(Nodes_ele,Nele)*s_0;
c_old = zeros(Num_Nodes,1);

%%%Evaluation of stress fibre concentration
while t<t_fin
	[gp_vec_x,w_gp_x]=Gauss_Points(n_gp);

gp_vec_y = gp_vec_x;
w_gp_y = w_gp_x;

n_gp_x = n_gp;
n_gp_y = n_gp;

N = zeros(n_gp_x,4,n_gp_y);
dN_dxi = zeros(n_gp_x,4,n_gp_y);
dN_deta = zeros(n_gp_x,4,n_gp_y);

N_bnd = zeros(2,4,2);
dN_dxi_bnd = zeros(2,4,2);
dN_deta_bnd = zeros(2,4,2);

N_mat = zeros(n_gp_x*Ndof,4*Ndof,n_gp_y);
N_bnd_mat = zeros(n_gp_x*Ndof,4*Ndof,n_gp_y);
coord_bnd = [-1 1];

%%Evaluation of bilinear shape functions
for j=1:n_gp_y
	eta = gp_vec_y(j);
	eta_bnd = coord_bnd(j);
	for i=1:n_gp_x
		xi = gp_vec_x(i);
		xi_bnd = coord_bnd(i);
		
		[N_vec,dN_dxi_vec,dN_deta_vec] = Shape_Function(xi,eta);
		
		[N_vec_bnd,dN_dxi_vec_bnd,dN_deta_vec_bnd] = Shape_Function(xi_bnd,eta_bnd);

		N(i,:,j) = N_vec(1,:,1); 
		dN_dxi(i,:,j) = dN_dxi_vec(1,:,1);
		dN_deta(i,:,j) = dN_deta_vec(1,:,1);

		N_bnd(i,:,j) = N_vec_bnd(1,:,1); 
		dN_dxi_bnd(i,:,j) = dN_dxi_vec_bnd(1,:,1);
		dN_deta_bnd(i,:,j) = dN_deta_vec_bnd(1,:,1);

		for k=1:Nodes_ele
			N_mat(i*2-1,k*2-1,j) = N(i,k,j);
			N_mat(i*2,k*2,j) = N(i,k,j);
			N_bnd_mat(i*2-1,k*2-1,j) = N_bnd(i,k,j);
			N_bnd_mat(i*2,k*2,j) = N_bnd(i,k,j);			
		end
		
	end
end

    % Cal = exp(-(t)/theta)*ones(Nele,1);          %Ad-hoc calcium signal, which initiates the system.

    %%%%%%%%%%%%%%%%%%%%%%%%% 
    %%Solve the first order differential equation to calculate concentration of stress fibres
    %%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(phi)
	eps_dot_old_inp_C = eps_dot_homo(:,ii);
	sigma_active_inp = sigma_active(:,ii);
	eta_old = eta_phi_old(:,ii);

	calcium_signal_solver;
	cal_out(:,ii) = Cal;
    cd RK45;
        main_solver_stress_fibre;         
    cd ..;

    eta_phi(:,ii) = eta_new(:);
end
% cal_out
% eta_phi

for  e = 1:Nele
	eta_homo_int(e) = trapz(phi,eta_phi(e,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = Elasticity_Matrix(E,nu);
condn = 1;
eps_r = 0.05;
del_u_max = 1;
iter = 0;

while condn > eps_r*del_u_max

	[Kt,K,F] = Stiffness_Matrix(N,dN_dxi,dN_deta,Nodes_ele,n_gp,Nodes,Elements,Ndof,Nele,Nodes_ele,...
							 Num_Nodes,C,N_bnd_mat,sigma_active,U_old_Eles,delta_t,kv_bar,epsdot_0,...
							 sig_max,sigma_ij,NI,eps_dot_homo_int,eta_homo_int,eps_homo_int,...
							 eps_old_homo_int,eps_0,eps_1,N_mat,xi_H_eles,F_rhs_old_inp_eles,lambda_s,...
							 shift_const,xi_0);

	H = K*U_old_inp - F;

	[Kt_red,H_red,U_Known,Disp_Flag,U_red_old] = Red_Stiffness_Matrix(Kt,H,Ndof,Num_Nodes,U_old_Eles,N_bnd_mat,bottom_nodes,top_nodes,right_nodes,left_nodes,U_old_inp);

	%%Solve for the displacements
	del_u = Kt_red\H_red;

	U_red = U_red_old - del_u;

	approx_err = abs((U_red - U_red_old)./U_red)*100;


	count_known = 1;
	count_unknown = 1;

	for i=1:Num_Nodes*Ndof
		if Disp_Flag(i) == 1
			U(i) = U_Known(count_known);
			count_known = count_known +1;
		else
			U(i) = U_red(count_unknown);
			count_unknown = count_unknown +1;
		end
	end

	condn = max(abs(U - U_old_inp));
	del_u_max = max(abs(U - U_old));

	U_old_inp = U;
	iter = iter + 1;

U_x = zeros(Num_Nodes,1);
U_y = zeros(Num_Nodes,1);

for i=1:Num_Nodes
	U_x(i) = U(i*2-1);
	U_y(i) = U(i*2);
end

%%%Split the Displacement vector element wise
U_Eles = zeros(Ndof*Nodes_ele,Nele);

for j=1:Nele
	ID = Elements(j,:);
	for i=1:Nodes_ele
		if Ndof == 1
			U_Eles(i,j) = U(ID(i));
		elseif Ndof == 2
			U_Eles(i*2-1,j) = U(ID(i)*2 -1);
			U_Eles(i*2,j) = U(ID(i)*2);
		end
	end
end

U_Dof_Eles = zeros(Ndof,Nodes_ele,Nele);

for e=1:Nele
	for j=1:Nodes_ele
		if Ndof == 1
			U_Dof_Eles(1,j,e) = U_Eles(j,e);
		elseif Ndof == 2
			U_Dof_Eles(1,j,e) = U_Eles(j*2-1,e);
			U_Dof_Eles(2,j,e) = U_Eles(j*2,e);
		end
	end
end

%%%Calculating new stretch
stretch_delta_new = zeros(Ndof*Num_Nodes,1);
for ii = 1:Ndof*Num_Nodes
	if abs(stretch_delta_old(ii)) < stretch_delta_max
		stretch_delta_new(ii) = stretch_delta_old(ii) + U(ii) - U_old(ii);
		if abs(stretch_delta_new(ii)) >= stretch_delta_max
			stretch_delta_new(ii) = stretch_delta_max*sign(stretch_delta_new(ii));
		end
	else
		stretch_delta_new(ii) = stretch_delta_max*sign(stretch_delta_new(ii));
	end
end
	
%%%Calculating new force applied
F_rhs_new = zeros(Num_Nodes*Ndof,1);
for ii = 1:Num_Nodes*Ndof
	F_rhs_new(ii) = lambda_s*stretch_delta_new(ii);
end

F_rhs_new_x = zeros(Num_Nodes,1);
F_rhs_new_y = zeros(Num_Nodes,1);

for i=1:Num_Nodes
	F_rhs_new_x(i) = F_rhs_new(i*2-1);
	F_rhs_new_y(i) = F_rhs_new(i*2);
end

%%%Split the Force vector element wise
F_rhs_eles = zeros(Ndof*Nodes_ele,Nele);

for j=1:Nele
	ID = Elements(j,:);
	for i=1:Nodes_ele
		if Ndof == 1
			F_rhs_eles(i,j) = F_rhs_new(ID(i));
		elseif Ndof == 2
			F_rhs_eles(i*2-1,j) = F_rhs_new(ID(i)*2 -1);
			F_rhs_eles(i*2,j) = F_rhs_new(ID(i)*2);
		end
	end
end

F_rhs_old_inp_eles = F_rhs_eles;

%%%Post processing
if Ndof == 1
	Eps_GP = zeros(2,n_gp,Nele);
elseif Ndof ==2
	Eps_GP = zeros(3,n_gp,Nele);
end
	
%%Stress and strain calculation 
for e=1:Nele
	ID = Elements(e,:);
	x_vec_ele(:,1) = Nodes(ID,1);
	y_vec_ele(:,1) = Nodes(ID,2);

	ID_Dof = Calc_ID(Ndof,ID);
	count_gp = 1;
	for j=1:n_gp
		for i=1:n_gp
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
			Eps_GP(:,count_gp,e) = B*U_Eles(:,e);
			Sigma_GP(:,count_gp,e) = C*Eps_GP(:,count_gp,e);
			count_gp = count_gp +1;
		end
	end
end

Eps_avg = zeros(3,Nele);
Sigma_avg = zeros(3,Nele);

for e = 1:Nele
	Eps_avg(1,e) = mean(Eps_GP(1,:,e));
	Eps_avg(2,e) = mean(Eps_GP(2,:,e));
	Eps_avg(3,e) = mean(Eps_GP(3,:,e));
	Sigma_avg(1,e) = mean(Sigma_GP(1,:,e));
	Sigma_avg(2,e) = mean(Sigma_GP(2,:,e));
	Sigma_avg(3,e) = mean(Sigma_GP(3,:,e));
end

Sigma_avg_11 = Sigma_avg(1,:);

%%calculate strain rate here
for ii=1:Nele
	Eps_dot(:,ii) = (Eps_avg(:,ii) - Eps_avg_old(:,ii))/delta_t;
end

%%Then calculate homogenised strain rate
for i=1:length(phi)
	eps_dot_homo(:,i) = Eps_dot(1,:)'*cos(phi(i))^2 + Eps_dot(2,:)'*sin(phi(i))^2 + Eps_dot(3,:)'*sin(2*phi(i));
	eps_homo(:,i) = Eps_avg(1,:)'*cos(phi(i))^2 + Eps_avg(2,:)'*sin(phi(i))^2 + Eps_avg(3,:)'*sin(2*phi(i));
end 

for i = 1:length(phi)
	for j = 1:Nele
			if eps_homo(j,i) <= 0
				f_eps = exp(-(eps_homo(j,i)/eps_0)^2);
			else
				f_eps = exp(-(eps_homo(j,i)/eps_0)^2) + (eps_homo(j,i)/eps_1)^2;
			end
            fac = 1/( 1 + ((shift_const)/( sqrt(shift_const^2 + 1)) ));
            sigma_active(j,i) = f_eps*fac*(eta_phi(j,i)*sig_max +   (eta_phi(j,i)*sig_max*(kv_bar*eps_dot_homo(j,i)/epsdot_0 + shift_const)/(sqrt( ((kv_bar)*(eps_dot_homo(j,i)/epsdot_0)  + shift_const)^2 + 1)) ));
    end
end

for e = 1:Nele

	sigma_ij(1,e) = (1/pi)*trapz(phi,sigma_active(e,:).*cos(phi).^2);
	sigma_ij(2,e) = (1/pi)*trapz(phi,sigma_active(e,:).*sin(phi).^2);
	sigma_ij(3,e) = (1/pi)*trapz(phi,(1/2)*sigma_active(e,:).*sin(2*phi));	

	eps_dot_homo_int(e) = (1/pi)*trapz(phi,eps_dot_homo(e,:));
	eps_homo_int(e) = (1/pi)*trapz(phi,eps_homo(e,:));

end

%%%Updating the concentration of low affinity integrins
main_solver_diffusion;

for iiii = 1:length(xi_H_new)
	if xi_H_new(iiii) < 0
		warning('bad bad bad')
	end
end

for iii = 1:length(xi_H_new)
	if xi_H_new(iii)/xi_0 >1
		warning('greater than 1')
	end
end


xi_H_dot_new = (xi_H_new - xi_H_old)/delta_t;

%%%Split the Displacement vector element wise
xi_H_eles = zeros(Nodes_ele,Nele);
xi_H_dot_new_eles = zeros(Nodes_ele,Nele);


for j=1:Nele
	ID = Elements(j,:);
	for i=1:Nodes_ele
		xi_H_eles(i,j) = xi_H_new(ID(i));
		xi_H_dot_new_eles(i,j) = xi_H_dot_new(ID(i));
	end
end

end
approx_err;
Eps_avg_old = Eps_avg;
eta_phi_old = eta_phi;
U_old_Eles = U_Eles;
U_old = U;
eps_old_homo = eps_homo;
eps_old_homo_int = eps_homo_int;
stretch_delta_old = stretch_delta_new;
xi_H_old = xi_H_new;

t_vec(index_t) = t;

if (t==0 || t==100 || t==500 || t==1000 || t==2000 || t==4000 || t==6000 || t==8000 || t==10000 )
	vtk_filename = strcat('res_',num2str(t),'.vtk');

	for e = 1:Nele
		eta_int(e) = trapz(phi,eta_phi(e,:)/pi);
	end

	Create_VTK(U,Num_Nodes,Nodes,Nele,Elements,Sigma_avg_11,eta_int,xi_H_new/xi_0,vtk_filename);
end

approx_err_plot(index_t) = norm(approx_err);


for ii = 1:length(phi)
	Cal_plot(index_t,ii) = mean(cal_out(:,ii));
end


for ii=1:length(phi)
	eta_plot(index_t,ii) = mean(eta_phi(:,ii));
end

for ii=1:length(phi)
	eps_plot(index_t,ii) = mean(eps_homo(:,ii));
end

for ii=1:length(phi)
	eps_dot_plot(index_t,ii) = mean(eps_dot_homo(:,ii));
end

for ii=1:length(phi)
	sigma_active_plot(index_t,ii) = mean(sigma_active(:,ii));
end


t = t+ delta_t
index_t = index_t+1;
end

% eta_plot_write = [t_vec' eta_plot];
% filename_eta = 'eta_with_FA';
% dlmwrite(filename_eta, eta_plot_write, 'delimiter', '\t', 'precision', '%f');

% eps_plot_write = [t_vec' eps_plot];
% filename_eps = 'eps_with_FA';
% dlmwrite(filename_eps, eps_plot_write, 'delimiter', '\t', 'precision', '%f');


% C_plot_write = [t_vec' Cal_plot];
% filename_eps = 'cal_with_FA';
% dlmwrite(filename_eps, C_plot_write, 'delimiter', '\t', 'precision', '%f');


% eta_tip = zeros(length(phi),1);
% for i = 1:length(phi)
% 	eta_tip(i) = eta_plot(r,i);
% end

% zer_vec = zeros(length(phi),1);

% data_polar = [ zer_vec zer_vec phi' eta_tip];
% filename_data = 'data_phi_1';
% dlmwrite(filename_data, data_polar, 'delimiter', '\t');


% % figure(1)
for ii = 6:length(phi)
	plot(t_vec,eta_plot(:,ii))
	hold all
end

% figure(2)
% for ii = 1:length(phi)
% 	plot(t_vec,eps_plot(:,ii))
% 	hold all
% end

% figure(3)
% for ii = 1:length(phi)
% 	plot(t_vec,eps_dot_plot(:,ii))
% 	hold all
% end

% figure(4)
% for ii = 1:length(phi)
% 	plot(t_vec,sigma_active_plot(:,ii))
% 	hold all
% end

% figure(5)
% for ii = 1:length(phi)
% 	plot(t_vec,Cal_plot(:,ii))
% 	hold all
% end

% figure(3)
% plot(t_vec,approx_err_plot);
