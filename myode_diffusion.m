%This function calculates the ODE needed by ode45 for integration.
function [xi_l_new,xi_alpha_nodes,F_cell_nodes] = myode_diffusion(stretch_delta_new,...
                                                                  lambda_s,Num_Nodes,Ndof,kT,xi_0,delta_mu)                                

stretch_delta_new_diff = stretch_delta_new;

F_cell_nodes = zeros(Num_Nodes,1);      %%%Force exerted by the cell on ECM, per node.
phi_ele_nodes = zeros(Num_Nodes,1);     %%%Energy calculated at each node.  
xi_alpha_nodes = zeros(Num_Nodes,1);    %%%Alpha calculated at each node
xi_alpha_dot_nodes = zeros(Num_Nodes*Ndof,1);

%calculating force exerted by the cell at every node
for i=1:Num_Nodes
    stretch_delta_new_eff = sqrt(stretch_delta_new_diff(i*2-1)^2 + stretch_delta_new_diff(i*2)^2);
    phi_ele_nodes = 1/2*lambda_s*stretch_delta_new_eff^2;
    stretch_delta_vec = [stretch_delta_new_diff(i*2-1);stretch_delta_new_diff(i*2)];

    xi_alpha_nodes(i) = exp(  (-delta_mu - phi_ele_nodes + lambda_s*stretch_delta_vec'*stretch_delta_vec)/kT  );
end

xi_l_new = xi_0./(1 + xi_alpha_nodes);