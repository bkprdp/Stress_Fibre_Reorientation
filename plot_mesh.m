function plot_mesh(Nodes,Elements,Nele,Nodes_ele)

	for e=1:Nele
		XX = 0;
		YY = 0;
		for i=1:Nodes_ele
			XX(i) = Nodes(Elements(e,i),1);
			YY(i) = Nodes(Elements(e,i),2);
			plot(XX,YY)
			hold on

			%Plot node numbers
			text(XX(i),YY(i),sprintf('%0.5g',Elements(e,i)));
		end
	end
