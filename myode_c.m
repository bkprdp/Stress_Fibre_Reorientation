%%%%Function to generate the differential equation for the growth of Calcium ions. t_c and y are the variables.

function diff_eqn = myode_c(t_C,y,lambda_f,lambda_b,s_inp,s_0)

diff_eqn = lambda_f*(s_inp/s_0).*(1-y) - lambda_b*y;