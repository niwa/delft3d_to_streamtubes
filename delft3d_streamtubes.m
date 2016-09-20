function [Nodes, Tubes] = delft3d_streamtubes(MdfFName, StreamtubesFName, ks, NoHorizTubes, NoVertTubes, CellsPerXs)
%DELFT3D_STREAMTUBES convert delft3D results to streamtubes
%   Detailed explanation goes here


%############### TO DO ################
% - Set filenames if not provided
% - Move flow stuff out so this function only processes a single flow - this
%   will make it much more generally applicable.
% - Pass H rather than MdfName


%% Set options if not provided
if ~exist('NoHorizTubes','var')
    NoHorizTubes = 10;
end
if ~exist('NoVertTubes','var')
    NoVertTubes = 5;
end
if ~exist('CellsPerXs','var')
    CellsPerXs = 2;
end

%% Read in model details
% model setup file (*.mdf)
MDF = delft3d_io_mdf('read' , MdfFName);
[ModelPath,MdfFName,~] = fileparts(MdfFName);

% model grid
Grid = delft3d_io_grd('read',fullfile(ModelPath,MDF.keywords.filcco));

% flowTS
BdyTS = bct_io('READ', fullfile(ModelPath, MDF.keywords.filbct));

% results file
H = vs_use(fullfile(ModelPath,['trim-',MdfFName,'.dat']));
T = vs_time(H);

%% Identify times at the end of a stationary flow period and associated flow
NFlows = 0;
for ii = 2:size(BdyTS.Table.Data,1)
    if BdyTS.Table.Data(ii,2) == BdyTS.Table.Data(ii-1,2)
        NFlows = NFlows+1;
        StationaryTimes(NFlows,1) = BdyTS.Table.Data(ii,1);
        StationaryFlows(NFlows,1) = BdyTS.Table.Data(ii,2);
    end
end

% Identify output times matching ends of stationary flow periods
NWanted = 0;
for ii = 1:size(T.t,1)
    if sum(StationaryTimes < T.t(ii)/T.tunit+1 & StationaryTimes > T.t(ii)/T.tunit-1) == 1
        NWanted = NWanted + 1;
        WantedTimes(NWanted,1) = ii;
        WantedFlows(NWanted,1) = StationaryFlows(StationaryTimes < T.t(ii)/T.tunit+1 & StationaryTimes > T.t(ii)/T.tunit-1);
    end
end

%% Read model results
% Get depth at specified times
S0 = vs_get(H,'map-const','DPS0',{2:Grid.nmax-1,2:Grid.mmax-1},'quiet');
WL = vs_get(H,'map-series',{WantedTimes},'S1',{2:Grid.nmax-1,2:Grid.mmax-1},'quiet');
Depth = cellfun(@(v) (v+S0),WL,'uniformoutput',false);
for ii = 1:length(Depth);
    Depth{ii}(isnan(Grid.cen.y)) = NaN;
    Depth{ii}(Depth{ii}<MDF.keywords.dryflc) = 0;
end;
clear S0 WL

% Get depth averaged velocity
U1 = vs_get(H,'map-series',{WantedTimes},'U1',{2:Grid.nmax-1,1:Grid.mmax-1,1},'quiet');
V1 = vs_get(H,'map-series',{WantedTimes},'V1',{1:Grid.nmax-1,2:Grid.mmax-1,1},'quiet');

% Checking I'm importing the right areas of the grid!
%KFU = vs_get(H,'map-series',{WantedTimes},'KFU',{2:Grid.nmax-1,1:Grid.mmax-1},'quiet');
%KFV = vs_get(H,'map-series',{WantedTimes},'KFV',{1:Grid.nmax-1,2:Grid.mmax-1},'quiet');

%% Identify dominant flow direction and prepare grids for processing
if abs(mean(mean(U1{end}))) > abs(mean(mean(V1{end})))
    % M direction dominant (use fliplr so processing is from Left to Right bank)
    Vel = cellfun(@fliplr,U1,'UniformOutput',false);
    Depth = cellfun(@fliplr,Depth,'UniformOutput',false);
    Xcor = fliplr(Grid.cor.x);
    Ycor = fliplr(Grid.cor.y);
else
    % N direction dominant (so transpose, fliplr not necessary)
    Vel = cellfun(@transpose,V1,'UniformOutput',false);
    Depth = cellfun(@transpose,Depth,'UniformOutput',false);
    Xcor = Grid.cor.x';
    Ycor = Grid.cor.y';
end

if mean(mean(Vel{end})) < 0
    % flow in negative direction so flip to process from upstream to downstream
    Vel = cellfun(@(v) -flipud(v),Vel,'UniformOutput',false);
    Depth = cellfun(@flipud,Depth,'UniformOutput',false);
    Xcor = flipud(Xcor);
    Ycor = flipud(Ycor);
end

% convert depths to cell faces
for FlowNo = 1:NWanted
    Depth{FlowNo} = [Depth{FlowNo}(:,1),Depth{FlowNo},Depth{FlowNo}(:,end)];
    Depth{FlowNo} = (Depth{FlowNo}(:,1:end-1) + Depth{FlowNo}(:,2:end)) /2;
end

% Tidy up
clear U1 V1


%% Build streamtubes
SelectedXs = [1:CellsPerXs:size(Depth{1},1),size(Depth{1},1)+1];
NoOfXs = size(SelectedXs,2);

Nodes = cell(NWanted,NoOfXs);
Tubes = cell(NWanted,NoOfXs);

% Loop through flows
for FlowNo = 1:NWanted
    % Loop through cross-sections
    XsCount = 0;
    for XsNo = SelectedXs
        XsCount = XsCount+1;
        XS_X     = Xcor(:,XsNo);
        XS_Y     = Ycor(:,XsNo);
        XS_Vel   = Vel{FlowNo}(:,XsNo);
        XS_Depth = Depth{FlowNo}(:,XsNo);
        [Nodes{FlowNo,XsCount},Tubes{FlowNo,XsCount}] = ...
            streamtubeXS(XS_X, XS_Y, XS_Vel, XS_Depth, ks, NoHorizTubes, NoVertTubes);
    end
end

%% Write out streamtubes file

% Loop through flows
for FlowNo = 1:NWanted
    Flow = WantedFlows(FlowNo);
    writeStreamtubes(Nodes, Tubes, StreamtubesFName, Flow, ks)
end

end

