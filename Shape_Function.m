function [N_vec,dN_dxi_vec,dN_deta_vec] = Shape_Function(xi,eta)

%This function will evaluate element level shape function
%and their derivatives.

% %number of gauss points
% n_gp_x = n_gp;
% n_gp_y = n_gp;

% N = zeros(n_gp_x,4,n_gp_y);
% dN_dxi = zeros(n_gp_x,4,n_gp_y);
% dN_deta = zeros(n_gp_x,4,n_gp_y);

%%Evaluation of bilinear shape functions
% for j=1:n_gp_y
% 	eta = gp_vec_y(j);
% 	for i=1:n_gp_x
% 		xi = gp_vec_x(i);
		N_vec = (1/4)*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];
		dN_dxi_vec = (1/4)*[-(1-eta) (1-eta) (1+eta) -(1+eta)];
		dN_deta_vec = (1/4)*[-(1-xi) -(1+xi) (1+xi) (1-xi)];

% 		N(i,:,j) = N_vec(1,:,1); 
% 		dN_dxi(i,:,j) = dN_dxi_vec(1,:,1);
% 		dN_deta(i,:,j) = dN_deta_vec(1,:,1);
% 	end
% end


