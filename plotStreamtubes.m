function plotStreamtubes(Nodes, Tubes)
%PLOTSTREAMTUBES plot streamtubes for a single flow
%   plotStreamtubes(Nodes, Tubes) plots streamtubes defined by Nodes and
%   Tubes.
%
%   Richard Measures 2016
%   
%   See also delft3d_streamtubes streamtubeXS

% Plot cross-sections
figure
for XsNo = 1:size(Tubes,2)
    for TubeNo = 1:(size(Tubes{1,XsNo},1)*size(Tubes{1,XsNo},2))
        plot3(Nodes{1,XsNo}(Tubes{1,XsNo}{TubeNo},1),...
              Nodes{1,XsNo}(Tubes{1,XsNo}{TubeNo},2),...
              Nodes{1,XsNo}(Tubes{1,XsNo}{TubeNo},3),'b-')
        hold on
    end
end

% Plot streamlines connecting cross-sections
NodeCoords = Nodes(1,:);
NodeCoords = permute(cell2mat(permute(NodeCoords,[1,3,2])),[3,1,2]);
for NodeNo = 1:size(NodeCoords,2)
    plot3(NodeCoords(:,NodeNo,1),NodeCoords(:,NodeNo,2),NodeCoords(:,NodeNo,3),'b-')
end
hold off
clear NodeCoords

end

