% @Author: pradeep
% @Date:   2015-12-09 11:33:56
% @Last Modified by:   pradeep
% @Last Modified time: 2015-12-09 11:37:27
function C = Elasticity_Matrix(E,nu)

	%Plane stress case

	C = E/(1-nu^2)*[1    nu     0         ;
					nu   1      0         ;
					0    0    (1-nu)/2   ];



	%Plane strain case 

	% C = E/((1+nu)(1-2*nu))*[ 1-nu     nu         0 			;
	%                          nu       1-nu       0 			;
	%                          0        0       (1-2*nu)/2   ];
	% 	