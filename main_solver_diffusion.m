%using xi_l at the current time step to calculate alpha, in order to update xi_H.
[xi_l_new,xi_alpha_nodes,F_cell_nodes] = myode_diffusion(stretch_delta_new,...
                                                         lambda_s,Num_Nodes,Ndof,kT,xi_0,delta_mu);                                
    
%calculate new xi_H
xi_H_new = xi_0 - xi_l_new;