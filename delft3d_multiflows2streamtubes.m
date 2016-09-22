function [Nodes, Tubes] = delft3d_multiflows2streamtubes(MdfFName, ...
    StreamtubesFName, ks, NoHorizTubes, NoVertTubes, CellsPerXs)
%DELFT3D_MULTIFLOWS2STREAMTUBES delft3d_streamtubes for a ramping flows
%   [Nodes, Tubes] = delft3d_multiflows2streamtubes(MdfFName, ...
%       StreamtubesFName, ks, NoHorizTubes, NoVertTubes, CellsPerXs)
%   
%   This function applies delft3d_streamtubes to every stable flow in a
%   ramping (stepped) flow timeseries e.g.
%   
%        |                *****x
%        |               *
%        |         *****x
%   Flow |        * 
%        |  *****x
%        |_______________________
%                  Time 
%   
%   where * marks flow time series and x marks output times used for
%   streamtube generation
%
%   Inputs:
%       MdfFname = String specifying filename (inc path if required) of 
%         delft3d model definition file (*.mdf). Note it is assumed model 
%         results are named consistently with this file and located in the 
%         same folder. If not supplied a user dialog box allows manual 
%         input.
%       StreamtubesFName = String specifying filename of output streamtubes
%         file. If not supplied a user dialog box allows manual input. Note
%         flow information is apended to the filename to distunguish
%         between the different flows in the range.
%       ks = roughness height in m. This must be manually input rather than
%         read from the model at the moment. Has a default value of 0.1 if
%         not supplied.
%       NoHorizTubes = Integer value specifying the number of tubes across 
%         the width of the river. Default value = 10.
%       NoVertTubes = Integer value specifying the number of layers of 
%         tubes in each vertical. Default value = 5.
%       CellsPerXs = Integer value specifying the number of delft3D cells
%         between each streamtubes cross-section. Default value = 2.
%
%   Outputs:
%       Nodes = cell array with size = [No_of_cross-sections,1] where each
%         cell is a matrix with size = [No_of_nodes,3] representing stream 
%         tube nodes at a single cross-section. Each row of the matrix 
%         represents 1 node, with the 3 columns specifying X, Y and Z
%         coordinates of the node (Z = depth below water surface).
%       Tubes = cell array with size = [No_of_cross-sections,1] where each
%         cell is a cell array with size = [NoVertTubes,NoHorizTubes]
%         representing the nodes associated with each tube at each
%         cross-section.
%
%   Richard Measures 2016
%
%   See also delft3d_streamtubes

%% Get File names if not supplied
if ~exist('MdfFName','var')
    [MdfFName,FilePath] = uigetfile('*.mdf','Select the delft3d model definition file');
    if isequal(MdfFName,0)
        error('User selected Cancel')
    end
    MdfFName = fullfile(FilePath,MdfFName);
    clear FilePath
end

if ~exist('StreamtubesFName','var')
    [StreamtubesFName,FilePath] = uiputfile('*.txt','Select file for streamtubes output', 'Streamtubes.txt');
    if isequal(StreamtubesFName,0)
        error('User selected Cancel')
    end
    StreamtubesFName = fullfile(FilePath,StreamtubesFName);
    clear FilePath
end

%% Set defaults if inputs not provided
if ~exist('NoHorizTubes','var')
    NoHorizTubes = 10;
end
if ~exist('NoVertTubes','var')
    NoVertTubes = 5;
end
if ~exist('CellsPerXs','var')
    CellsPerXs = 2;
end
if ~exist('ks','var') % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    warning('Setting ks to default value of 0.1m')
    ks = 0.1;
end

%% Read in model details
% model setup file (*.mdf)
MDF = delft3d_io_mdf('read' , MdfFName);
[ModelPath,ModelName,~] = fileparts(MdfFName);

% flowTS
BdyTS = bct_io('READ', fullfile(ModelPath, MDF.keywords.filbct));

% results file
H = vs_use(fullfile(ModelPath,['trim-',ModelName,'.dat']));
T = vs_time(H);

%% Identify output times matching ends of stationary flow periods

% Identify the timeseries inflow boundary condition
NoOfTotDisBdy = 0;
for ii = 1:BdyTS.NTables
    if strcmp(BdyTS.Table(ii).Parameter(2).Name(1:15),'total discharge');
        NoOfTotDisBdy = NoOfTotDisBdy + 1;
        InBdyNo = ii;
    end
end
if NoOfTotDisBdy ~= 1
    error('Too many total discharge boundaries for program to cope!')
end

% Identify times at the end of a stationary flow period and associated flow
NFlows = 0;
for ii = 2:size(BdyTS.Table(InBdyNo).Data,1)
    if BdyTS.Table(InBdyNo).Data(ii,2) == BdyTS.Table(InBdyNo).Data(ii-1,2)
        NFlows = NFlows+1;
        StationaryTimes(NFlows,1) = BdyTS.Table(InBdyNo).Data(ii,1);
        StationaryFlows(NFlows,1) = BdyTS.Table(InBdyNo).Data(ii,2);
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

%% Loop over wanted times and generate streamtubes for each one

[OutPath,OutName,OutExt] = fileparts(StreamtubesFName);
Nodes = cell(NWanted,1);
Tubes = cell(NWanted,1);

% Loop through flows
for FlowNo = 1:NWanted
    Flow = WantedFlows(FlowNo);
    OutFullFile = strrep(sprintf('%s_Flow=%07.3f',OutName, Flow),'.','-');
    OutFullFile = fullfile(OutPath,[OutFullFile,OutExt]);
    
    [Nodes{FlowNo,1}, Tubes{FlowNo,1}] = ...
        delft3d_streamtubes(MdfFName, WantedTimes(FlowNo), OutFullFile, ...
                            ks, NoHorizTubes, NoVertTubes, CellsPerXs, ...
                            Flow);
end


end

