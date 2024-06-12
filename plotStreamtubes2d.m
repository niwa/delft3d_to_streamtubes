function FigureH = plotStreamtubes2d(Nodes, Tubes, Mask)
%PLOTSTREAMTUBES plot map of surface streamtubes for a single flow
%   plotStreamtubes2d(Nodes, Tubes) plots streamtubes defined by Nodes and
%   Tubes.
%
%   Richard Measures, NIWA, 2016
%   
%   See also plotStreamtubes delft3d_streamtubes streamtubeXS

FigureH = figure;

NoOfNodes = size(Nodes{1,1},1);
NoOfVerts = size(Tubes{1,1},2);
NoOfLayers = size(Tubes{1,1},1);

% Plot streamlines
NodeCoords = Nodes(:);
NodeCoords = permute(cell2mat(permute(NodeCoords,[2,3,1])),[3,1,2]);
for NodeNo = 1:NoOfNodes
    if NodeCoords(1,NodeNo,3) == 0 %this plots only the first layer
        plot(NodeCoords(:,NodeNo,1),NodeCoords(:,NodeNo,2),'k-');
        hold on
    end
end
clear NodeCoords

% Plot cross-sections

for XsNo = 1:size(Tubes,1)
    Tube=Tubes{XsNo,1};
    MaskTube = Mask{XsNo,1};
    LayerNo = 1; %use only the first 
    for VertNo = 1:NoOfVerts
        NodeNumber = Tube{LayerNo,VertNo}; % %Node number (4 values that define one tube
        MaskValue = MaskTube{LayerNo,VertNo};
        
        if(MaskValue==1)
            plot(Nodes{XsNo,1}(NodeNumber,1),...
                 Nodes{XsNo,1}(NodeNumber,2),...
                 'b-');
        else
            plot(Nodes{XsNo,1}(NodeNumber,1),...
                 Nodes{XsNo,1}(NodeNumber,2),...
                 'r-');
        end
    end
    hold on
end

hold off

% adjust formatting
axis equal
xlabel('Easting [m]')
ylabel('Northing [m]')

end

