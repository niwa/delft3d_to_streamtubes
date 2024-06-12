function writeStreamtubes(Nodes, Tubes, OutFName, Flow, ks, Mask, OutFNameMask)
%WRITESTREAMTUBES output streamtubes to text file
%   writeStreamtubes(Nodes, Tubes, OutFName)
%
%   Richard Measures, Gu Stecca, NIWA
%
%   See also delft3d_streamtubes

%% Extract some basic parameters
NoOfXs     = size(Nodes,1);
NoOfVerts  = size(Tubes{1,1},2);
NoOfLayers = size(Tubes{1,1},1);
NoOfNodes  = size(Nodes{1,1},1);

%% Flexibility to account for either single or variable ks
if length(ks) == 1
    ks = repmat(ks,NoOfVerts+1,1);
end

%% Open file for write access
FID=fopen(OutFName,'w');

%% Write header info
fprintf(FID, 'Drift Model Output\r\n');
fprintf(FID, 'number of x-secs: %i\r\n', NoOfXs);
fprintf(FID, 'number of flow tubes: %i\r\n', NoOfVerts);
fprintf(FID, 'number of VertTubes: %i\r\n', NoOfLayers);

%% Write cross-section data
for XsNo = 1:NoOfXs
    % Cross-section header
    fprintf(FID, '%-7i %1.5f\r\n',XsNo,Flow);
    fprintf(FID, '%-7.5f ',ks);
    fprintf(FID, '\r\n');
    % Nodes
    for NodeNo = 1:NoOfNodes
        fprintf(FID, '%-7i %-7.2f %-7.2f % -7.6f\r\n', NodeNo, Nodes{XsNo,1}(NodeNo,:));
    end
    % Tubes
    for VertNo = 1:NoOfVerts
        for LayerNo = 1:NoOfLayers
            fprintf(FID, '%-7i ',Tubes{XsNo,1}{LayerNo,VertNo});
            fprintf(FID, '\r\n');
        end
    end
end

%% Close file
fclose(FID);


if(exist(Mask))
    
    FIDM=fopen(OutFNameMask,'w');
    
    %% Write header info
    fprintf(FIDM, 'Drift Model Output\r\n');
    fprintf(FIDM, 'number of x-secs: %i\r\n', NoOfXs);
    fprintf(FIDM, 'number of flow tubes: %i\r\n', NoOfVerts);
    fprintf(FIDM, 'number of VertTubes: %i\r\n', NoOfLayers);
    
    %% Write cross-section data
    for XsNo = 1:NoOfXs
        % Cross-section header
        fprintf(FIDM, '%-7i %1.5f\r\n',XsNo,Flow);
        % Masks
        for VertNo = 1:NoOfVerts
            for LayerNo = 1:NoOfLayers
                fprintf(FIDM, '%-7i ',Mask{XsNo,1}{LayerNo,VertNo});
                fprintf(FIDM, '\r\n');
            end
        end
    end
    fclose(FIDM);
    
end


end

