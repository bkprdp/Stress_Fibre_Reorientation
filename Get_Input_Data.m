function [  Nodes,Elements,Left_Edge,Right_Edge,Top_Edge,Bottom_Edge,Node_BR,Node_BL,Node_TR,Node_TL] = Get_Input_Data(filename)

fid = fopen (filename, 'r');

line = fgetl (fid);
while (all (line ~= -1) || isempty (line))
  % read the Nodes
  if (strcmpi (line, '<Nodes>'))
      line = fgetl (fid);
      Nodes = [];
      EndNodes = 0;
      while EndNodes == 0
        NewNode= strread(line,'%f ','delimiter', '; ')';
        Nodes = [Nodes; NewNode(2:end)];
        line = fgetl (fid);
        if (strcmpi (line, '</Nodes>'))
          EndNodes = 1;
        end
      end
  end

  % read the elements
  if (strcmpi (line, '<Elements>'))
    line = fgetl (fid);
    Elements = [];
    EndElements = 0;
    while EndElements == 0
        NewElement= strread(line,'%d ','delimiter', '; ')';
        Elements = [Elements; NewElement(2:end)+1];
        line = fgetl (fid);
        if (strcmpi (line, '</Elements>'))
          EndElements = 1;
        end
      end
  end


  %Read nodes on left edge
  if (strcmpi (line, '<NodeGroup name="Left">'))
    line = fgetl (fid) (2:end-1);
    Left_Edge = textscan(line,'%d'){1}+1;
  end


  % Read nodes at Right edge
  if (strcmpi (line, '<NodeGroup name="Right">'))
    line = fgetl (fid) (2:end-1);
    Right_Edge = textscan(line,'%d'){1}+1;
  end

  % Read nodes on top edge
  if (strcmpi (line, '<NodeGroup name="Top">'))
    line = fgetl (fid) (2:end-1);
    Top_Edge = textscan(line,'%d'){1}+1;
  end

  % Read nodes at bottom edge
  if (strcmpi (line, '<NodeGroup name="Bottom">'))
    line = fgetl (fid) (2:end-1);
    Bottom_Edge= textscan(line,'%d'){1}+1;
  end

  %Read Corner nodes
  if (strcmpi (line, '<Node name="Corner_BR">'))
    line = fgetl (fid) (2:end-1);
    Node_BR= textscan(line,'%d'){1}+1;
  end

  if (strcmpi (line, '<Node name="Corner_BL">'))
    line = fgetl (fid) (2:end-1);
    Node_BL= textscan(line,'%d'){1}+1;
  end

  if (strcmpi (line, '<Node name="Corner_TR">'))
    line = fgetl (fid) (2:end-1);
    Node_TR= textscan(line,'%d'){1}+1;
  end
  
  if (strcmpi (line, '<Node name="Corner_TL">'))
    line = fgetl (fid) (2:end-1);
    Node_TL= textscan(line,'%d'){1}+1;
  end
  line = fgetl (fid);
end

fclose (fid);

end
