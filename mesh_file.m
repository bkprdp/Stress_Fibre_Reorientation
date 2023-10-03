function [Nodes,Elements,Nele,Nodes_ele,Num_Nodes,bottom_nodes,top_nodes,right_nodes,left_nodes] = mesh_file

            % Geometry
            l_x = 12.5e-6;
            l_y = 12.5e-6;
            div_x = 5;
            div_y = 5;
            
            nno=0;

            
            % creating nodes
            for i=0:div_y
                left_nodes(i+1) = nno+1;
                for j=0:div_x
                    nno = j+1+(div_x+1)*i;
                    x = l_x/div_x*j;
                    y = l_y/div_y*i;
                    Nodes(nno,:) = [x y]; 
                    if i==0 
                        bottom_nodes(j+1) = nno;
                    end
                    if i==div_y
                        top_nodes(j+1) = nno;
                    end
                end
                right_nodes(i+1) = nno;
            end
            
            
            % creating elements
            for i=0:div_y-1
                for j=0:div_x-1
                    eno = j+1+(div_x)*i;
                    n1 = j+1+(div_x+1)*i;
                    n2 = j+1+(div_x+1)*i + 1;
                    n4 = j+1+(div_x+1)*(i+1);
                    n3 = j+1+(div_x+1)*(i+1) + 1;
                    Elements(eno,:) = [n1 n2 n3 n4];
                end
            end

            Num_Nodes = nno;
            Nele = eno;
            Nodes_ele = 4;