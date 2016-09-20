function writeStreamtubes(Nodes, Tubes, OutFName, Flow, ks)
%WRITESTREAMTUBES output streamtubes to text file
%   writeStreamtubes(Nodes, Tubes, OutFName)
%
%   Richard Measures 2016
%
%   See also delft3d_streamtubes

%% Extract some basic parameters
NoOfXs     = size(Nodes,2);
NoOfVerts  = size(Tubes{1,1},2);
NoOfLayers = size(Tubes{1,1},1);

%% Flexibility to account for either single or variable ks
if length(ks) == 1
    ks = repmat(ks,NoOfVerts+1,1);
end

%% Open file for write access
FID=fopen(OutFName,'w');

%% Write header info
fprintf(FID, 'number of x-secs: %i\n', NoOfXs);
fprintf(FID, 'number of flow tubes: %i\n', NoOfVerts);
fprintf(FID, 'number of VertTubes: %i\n', NoOfLayers);

%% Write cross-section data
for XsNo = 1:NoOfXs
    % Cross-section header
    fprintf(FID, '%-7i %1.5f\n',XsNo,Flow);
    fprintf(FID, '%-7.5f ',ks);
    fprintf(FID, '\n');
    % Nodes
    for NodeNo = 1:size(Nodes{1,XsNo},1)
        fprintf(FID, '%-7i %-7.2f %-7.2f %-7.6f\n', NodeNo, Nodes{1,XsNo}(NodeNo,:));
    end
    % Tubes
    for VertNo = 1:NoOfVerts
        for LayerNo = 1:NoOfLayers
            fprintf(FID, '%-7i ',Tubes{1,XsNo}{LayerNo,VertNo});
            fprintf(FID, '\n');
        end
    end
end

%% Close file
fclose(FID);

end

