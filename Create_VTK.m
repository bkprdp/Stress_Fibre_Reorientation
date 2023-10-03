% @Author: pradeep
% @Date:   2015-12-09 15:50:23
% @Last Modified by:   pradeep
% @Last Modified time: 2015-12-16 17:11:10
function Create_VTK(U,Num_Nodes,Nodes,Nele,Elements,Sigma_avg,eta_new,xi_H,filename_VTK)

% Open the file for writing
fid = fopen(filename_VTK,'w');

% print header
fprintf(fid,'# vtk DataFile Version 3.1 \n');
fprintf(fid,'This is the vtk file for the elasticity problem \n');
fprintf(fid,'ASCII \n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
fprintf(fid,'POINTS %d FLOAT \n',Num_Nodes);

% print Nodal coordinates
for i=1:Num_Nodes
    fprintf(fid,'%.12g %.12g 0.\n',Nodes(i,1),Nodes(i,2));
end

% print elements
fprintf(fid,'CELLS %12d %12d \n',Nele,5*Nele);
for e=1:Nele
    fprintf(fid,'4 %12d %12d %12d %12d \n',Elements(e,1)-1,Elements(e,2)-1,Elements(e,3)-1,Elements(e,4)-1);
end

fprintf(fid,'CELL_TYPES %12d \n',Nele);
for e=1:Nele
    fprintf(fid,'9  ');
end
fprintf(fid,'\n');

fprintf(fid,'POINT_DATA %12d \n',Num_Nodes);
fprintf(fid,'SCALARS u_x FLOAT \n');
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:Num_Nodes
	fprintf(fid, '%e \n',U(i*2-1));
end

fprintf(fid,'SCALARS u_y FLOAT \n');
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:Num_Nodes
	fprintf(fid, '%e \n',U(i*2));
end

fprintf(fid,'SCALARS xi_H FLOAT \n');
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:Num_Nodes
	fprintf(fid, '%e \n',xi_H(i));
end


% fprintf(fid,'VECTORS xi_H FLOAT \n');
% for i=1:Num_Nodes
% 	fprintf(fid, '%e %e %e \n',xi_H(i*2-1)/sqrt(2),xi_H(i*2)/sqrt(2),0);
% end

fprintf(fid,'CELL_DATA %12d \n',Nele);
fprintf(fid,'SCALARS sigma_x FLOAT \n');
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:Nele
	fprintf(fid, '%e \n',Sigma_avg(1,i));
end

fprintf(fid,'SCALARS eta FLOAT \n');
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:Nele
	fprintf(fid, '%e \n',eta_new(i));
end
