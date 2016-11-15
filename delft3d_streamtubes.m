function [Nodes, Tubes] = delft3d_streamtubes(MdfFName, OutputTimeID, StreamtubesFName, ks, NoHorizTubes, NoVertTubes, CellsPerXs, TotalFlow)
%DELFT3D_STREAMTUBES convert delft3D results to streamtubes
%   [Nodes, Tubes] = ...
%       delft3d_streamtubes(MdfFName, OutputTimeID, StreamtubesFName, ...
%                           ks, NoHorizTubes, NoVertTubes, CellsPerXs)
%
%   Inputs:
%       MdfFname = String specifying filename (inc path if required) of 
%         delft3d model definition file (*.mdf). Note it is assumed model 
%         results are named consistently with this file and located in the 
%         same folder. If not supplied a user dialog box allows manual 
%         input.
%       OutputTimeID = Integer refering to the timestep of the trim file
%         from which results should be extracted. Default is the last
%         timestep in the trim file.
%       StreamtubesFName = String specifying filename of output streamtubes
%         file. If not supplied a user dialog box allows manual input.
%       ks = roughness height in m. This must be manually input rather than
%         read from the model at the moment. Has a default value of 0.1 if
%         not supplied.
%       NoHorizTubes = Integer value specifying the number of tubes across 
%         the width of the river. Default value = 20.
%       NoVertTubes = Integer value specifying the number of layers of 
%         tubes in each vertical. Default value = 5.
%       CellsPerXs = Integer value specifying the number of delft3D cells
%         between each streamtubes cross-section. Default value = 2.
%       TotalFlow = Optional user specifed total flow. If not supplied
%         calculated mean cross-section flow is used instead.
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
%   Richard Measures, NIWA, 2016
%
%   See also streamtubeXS delft3d_multiflows2streamtubes

%% ############### TO DO ################
% - Improve calculation of depth at cell faces taking into account dpsopt
%   and dpuopt... confusing!
% - Deal with spatially varying roughness and get roughness from model

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
    NoHorizTubes = 20;
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
[ModelPath,MdfFName,~] = fileparts(MdfFName);

% model grid
Grid = delft3d_io_grd('read',fullfile(ModelPath,MDF.keywords.filcco));

% results file
H = vs_use(fullfile(ModelPath,['trim-',MdfFName,'.dat']));
T = vs_time(H);

%% Read model results

% Set default output time if not provided
if ~exist('OutputTimeID','var')
    OutputTimeID = T.nt_loaded;
end

% Get depth at specified time
S0 = vs_get(H,'map-const','DPS0',{2:Grid.nmax-1,2:Grid.mmax-1},'quiet');
WL = vs_get(H,'map-series',{OutputTimeID},'S1',{2:Grid.nmax-1,2:Grid.mmax-1},'quiet');
Depth = S0 + WL;
Depth(isnan(Grid.cen.y)) = NaN;
Depth(Depth<MDF.keywords.dryflc) = 0;
clear S0 WL

% Get depth averaged velocity
U1 = vs_get(H,'map-series',{OutputTimeID},'U1',{2:Grid.nmax-1,1:Grid.mmax-1,1},'quiet');
V1 = vs_get(H,'map-series',{OutputTimeID},'V1',{1:Grid.nmax-1,2:Grid.mmax-1,1},'quiet');

% Checking I'm importing the right areas of the grid!
%KFU = vs_get(H,'map-series',{OutputTimeID},'KFU',{2:Grid.nmax-1,1:Grid.mmax-1},'quiet');
%KFV = vs_get(H,'map-series',{OutputTimeID},'KFV',{1:Grid.nmax-1,2:Grid.mmax-1},'quiet');

%% Identify dominant flow direction and prepare grids for processing

if abs(mean(mean(U1))) > abs(mean(mean(V1)))
    % M direction dominant (use flipud so processing is from Left to Right bank)
    Vel   = flipud(U1);
    Depth = flipud(Depth);
    Xcor  = flipud(Grid.cor.x);
    Ycor  = flipud(Grid.cor.y);
else
    % N direction dominant (so transpose, flipud not necessary)
    Vel   = V1';
    Depth = Depth';
    Xcor  = Grid.cor.x';
    Ycor  = Grid.cor.y';
end

if mean(mean(Vel)) < 0
    % flow in negative direction so fliplr to process from upstream to downstream
    Vel   = -flipud(Vel);
    Depth = flipud(Depth);
    Xcor  = flipud(Xcor);
    Ycor  = flipud(Ycor);
end

% convert depths to cell faces
Depth = [Depth(:,1),Depth,Depth(:,end)];
Depth = (Depth(:,1:end-1) + Depth(:,2:end)) /2;

% Tidy up
clear U1 V1

%% Build streamtubes

SelectedXs = [1:CellsPerXs:size(Depth,2)-1,size(Depth,2)];
NoOfXs = size(SelectedXs,1);

Nodes = cell(NoOfXs,1);
Tubes = cell(NoOfXs,1);
Flows = nan(NoOfXs,1);

% Loop through cross-sections
XsCount = 0;
for XsNo = SelectedXs
    XsCount = XsCount+1;
    XS_X     = Xcor(:,XsNo);
    XS_Y     = Ycor(:,XsNo);
    XS_Vel   = Vel(:,XsNo);
    XS_Depth = Depth(:,XsNo);
    [Nodes{XsCount,1},Tubes{XsCount,1},Flows(XsCount,1)] = ...
        streamtubeXS(XS_X, XS_Y, XS_Vel, XS_Depth, ks, ...
                     NoHorizTubes, NoVertTubes);
end
clear XS_X XS_Y XS_Vel XS_Depth XS_Count

%% Write out streamtubes file

MeanFlow = mean(Flows);
if ~exist('TotalFlow','var')
    TotalFlow = MeanFlow;
else
    if (TotalFlow > (MeanFlow * 1.1)) || (TotalFlow < (MeanFlow * 0.9))
        warning('TotalFlow specified (%.3f m^3/s) is significantly different to meanflow at cross sections (%.3f m^3/s)',TotalFlow,MeanFlow)
    end
end
writeStreamtubes(Nodes, Tubes, StreamtubesFName, TotalFlow, ks)

end

