function FigureH = plotStreamtubes(Nodes, Tubes)
%PLOTSTREAMTUBES plot streamtubes for a single flow
%   plotStreamtubes(Nodes, Tubes) plots streamtubes defined by Nodes and
%   Tubes.
%
%   Richard Measures, NIWA, 2016
%   
%   See also plotStreamtubes2d delft3d_streamtubes streamtubeXS

FigureH = figure;

% Plot streamlines
NodeCoords = Nodes(:);
NodeCoords = permute(cell2mat(permute(NodeCoords,[2,3,1])),[3,1,2]);
for NodeNo = 1:size(NodeCoords,2)
    plot3(NodeCoords(:,NodeNo,1),NodeCoords(:,NodeNo,2),NodeCoords(:,NodeNo,3),'b-')
    hold on
end
clear NodeCoords

% Plot cross-sections
for XsNo = 1:size(Tubes,1)
    for TubeNo = 1:(size(Tubes{XsNo,1},1)*size(Tubes{XsNo,1},2))
        patch(Nodes{XsNo,1}(Tubes{XsNo,1}{TubeNo},1),...
              Nodes{XsNo,1}(Tubes{XsNo,1}{TubeNo},2),...
              Nodes{XsNo,1}(Tubes{XsNo,1}{TubeNo},3),'w','FaceAlpha',0.8)
    end
end
hold off

% adjust formatting
daspect([20 20 1])
xlabel('Easting [m]')
ylabel('Northing [m]')
zlabel('Depth [m]')
whitebg('k')

end

