function [XS_Nodes,XS_Tubes,TotalFlow,XS_Mask,XS_ks_Node] = streamtubeXS(XS_X, XS_Y, XS_Vel, XS_Depth, XS_ks, NoHorizTubes, ...
                                                     NoVertTubes, MaxDryCellsInTube)
%streamtubeXS Generate streamtubes for a single cross-section
%   [Nodes,Tubes,TotalFlow] = streamtubeXS(XS_X, XS_Y, XS_Vel, XS_Depth,...
%                                          XS_ks, NoHorizTubes, NoVertTubes)
%
%   Richard Measures, Gu Stecca, NIWA
%
%   See also delft3d_streamtubes


%% First, some tidying of the input data

% Remove any negative or nan velocities
XS_Vel(XS_Vel<0|isnan(XS_Vel)) = 0;

% Remove any negative or nan depths
XS_Depth(XS_Depth<0|isnan(XS_Depth)) = 0;
DryCells = XS_Depth<=0;

%% Calcuate some basic cross-section properties

% Calculate cell width
CellWidth = sqrt((XS_X(1:end-1)-XS_X(2:end)).^2 + (XS_Y(1:end-1)-XS_Y(2:end)).^2);
CellWidth(isnan(CellWidth)) = 0;

% Calculate cell flow
CellFlow = XS_Vel .* XS_Depth .* CellWidth;
CumulativeFlow = [0;cumsum(CellFlow)];
TotalFlow = sum(CellFlow);

% Find Banks
FlowingCells = CellFlow > 0;
FlowingEdges = [false;FlowingCells]|[FlowingCells;false];

%% Divide XS horizontally into streamtubes

% 1st make sure CumulativeFlow increases so we can use it in interp1
CumFlowOfFlowingEdges = CumulativeFlow(FlowingEdges);
while sum(CumFlowOfFlowingEdges(2:end) == CumFlowOfFlowingEdges(1:end-1))>0
    CumFlowOfFlowingEdges([false;(CumFlowOfFlowingEdges(2:end) == CumFlowOfFlowingEdges(1:end-1))]) = ...
        CumFlowOfFlowingEdges([false;(CumFlowOfFlowingEdges(2:end) == CumFlowOfFlowingEdges(1:end-1))]) + 1e-8;
end

TubeEdgePos = interp1(CumFlowOfFlowingEdges,...
                      find(FlowingEdges),...
                      [0,(TotalFlow/NoHorizTubes)*(1:NoHorizTubes-1),max(CumFlowOfFlowingEdges)]');

TubeEdgeXYDVks = [interp1([XS_X, XS_Y], TubeEdgePos),...
                  interp1([XS_Depth(1)  , XS_Vel(1)  , XS_ks(1)  ;...
                          XS_Depth     , XS_Vel     , XS_ks     ;...
                          XS_Depth(end), XS_Vel(end), XS_ks(end)],...
                          (TubeEdgePos+0.5))];

TubeEdgePos_L= floor(TubeEdgePos(1:end-1));
TubeEdgePos_R= ceil(TubeEdgePos(2:end))-1;

MaskTubes = ones(size(TubeEdgeXYDVks, 1)-1, 1);
 
for TubeNo = 1:size(MaskTubes, 1)
    NDry = sum(DryCells(TubeEdgePos_L(TubeNo):TubeEdgePos_R(TubeNo)));
    MaskTubes(TubeNo) = 1*(NDry<=MaxDryCellsInTube);
end
 
%% Divide XS vertically into streamtubes

% Calculate vertical divisions
NoBins = 100; % hardcoded!
TubeEdgeZ = zeros(NoHorizTubes+1,NoVertTubes+1);
for EdgeNo = 2:NoHorizTubes
    zBins = [1e-6,(1:NoBins)*(TubeEdgeXYDVks(EdgeNo,3)/NoBins)]';
    IntegratedVel = log(11*zBins/TubeEdgeXYDVks(EdgeNo,5)).*zBins;
    TubeEdgeZ(EdgeNo,2:end-1) = interp1(IntegratedVel,zBins,(1:NoVertTubes-1)*IntegratedVel(end)/NoVertTubes);
    TubeEdgeZ(EdgeNo,1:end-1) = -TubeEdgeXYDVks(EdgeNo,3) + TubeEdgeZ(EdgeNo,1:end-1);
end

%% Convert to nodes and tubes (polygons)

% Nodes
% Node numbering = 1                                                  for Left Bank (EdgeNo=1)
%                = (EdgeNo-2)*(NoVertTubes+1) + (2:(NoVertTubes+2))   for intermediate verticals (EdgeNo = 2:NoHorizTubes)
%                = (NoHorizTubes-1) * (NoVertTubes+1) + 2             for Right Bank (EdgeNo = NoHorizTubes+1)
XS_Nodes = nan((NoHorizTubes-1) * (NoVertTubes+1) + 2, 3);
XS_Nodes(1,:) = [TubeEdgeXYDVks(1,[1,2]),0]; % Left bank
for EdgeNo = 2:NoHorizTubes;
    XS_Nodes((EdgeNo-2)*(NoVertTubes+1) + (2:(NoVertTubes+2)), :) = ...
        [repmat(TubeEdgeXYDVks(EdgeNo,[1,2]),NoVertTubes+1,1),TubeEdgeZ(EdgeNo,:)'];
end
XS_Nodes(end,:) = [TubeEdgeXYDVks(end,[1,2]),0]; % Right bank

% Tubes
XS_Tubes = cell(NoVertTubes,NoHorizTubes);
XS_Mask = cell(NoVertTubes,NoHorizTubes);
% Left Bank Tubes
for LayerNo = 1:NoVertTubes
    XS_Tubes{LayerNo,1} = [1,LayerNo+(1:2),1];
end
% Middle Tubes
for VertNo = 2:NoHorizTubes-1
    for LayerNo = 1:NoVertTubes
        XS_Tubes{LayerNo,VertNo} = [(VertNo-2)*(NoVertTubes+1)+LayerNo+1,...
                                    (VertNo-1)*(NoVertTubes+1)+LayerNo+[1,2],...
                                    (VertNo-2)*(NoVertTubes+1)+LayerNo+2];
    end
end
% Right Bank Tubes
for LayerNo = 1:NoVertTubes
    XS_Tubes{LayerNo,end} = [(NoHorizTubes-2)*(NoVertTubes+1)+LayerNo+1,...
                             (NoHorizTubes-1)*(NoVertTubes+1)+2,...
                             (NoHorizTubes-2)*(NoVertTubes+1)+LayerNo+[2,1]];
end

for VertNo = 1:NoHorizTubes
    for LayerNo = 1:NoVertTubes
        XS_Mask{LayerNo,VertNo} = MaskTubes(VertNo,1);
    end
end

% ks values in the nodes (only the lower layer can be interpreted)
XS_ks_Node = [XS_ks(TubeEdgePos_L); XS_ks(TubeEdgePos_L(end))];


%% Plot the XS (testing purposes only)

% % Plot vs Y
% plot(XS_Y,-[XS_Depth(1);(XS_Depth(1:end-1)+XS_Depth(2:end))/2;XS_Depth(end)],'r.-')
% hold on
% plot(XS_Nodes(:,2),XS_Nodes(:,3),'go')
% for TubeNo = 1:NoHorizTubes*NoVertTubes
%     plot(XS_Nodes(XS_Tubes{TubeNo},2),XS_Nodes(XS_Tubes{TubeNo},3),'b-')
% end
% hold off
% 
% % Plot vs X
% plot(XS_X,-[XS_Depth(1);(XS_Depth(1:end-1)+XS_Depth(2:end))/2;XS_Depth(end)],'r.-')
% hold on
% plot(XS_Nodes(:,1),XS_Nodes(:,3),'go')
% for TubeNo = 1:NoHorizTubes*NoVertTubes
%     plot(XS_Nodes(XS_Tubes{TubeNo},1),XS_Nodes(XS_Tubes{TubeNo},3),'b-')
% end
% hold off

end

