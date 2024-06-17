function [Nodes, Tubes, Mask] = delft3d_streamtubes(MdfFName, ...
    OutputTimeID, StreamtubesFName, MaskFName, ks, NoHorizTubes, ...
    NoVertTubes, CellsPerXs, TotalFlow, MaxDryCellsInTube)
%DELFT3D_STREAMTUBES convert delft3D results to streamtubes
%   [Nodes, Tubes, Mask] = ...
%       delft3d_streamtubes(MdfFName, OutputTimeID, StreamtubesFName, ...
%                           MaskFName, ks, NoHorizTubes, NoVertTubes, ...
%                             CellsPerXs, TotalFlow, MaxDryCellsInTube)
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
%       MaskFName = String specifying filename of output streamtubes mask
%         file. If not supplied a user dialog box allows manual input.
%       ks = roughness height in m. This can be specified manually. If no
%           value is specified (that is, [] must be entered as an input
%           value), then it is read from the model (check if it is set as a
%           spatially constant value in MDF.keywords.ccofu or as a
%           spatially varying set in the .rgh file). If none of them are
%           found (which in principle cannot be), then the default value
%           ks = 0.1.
%       NoHorizTubes = Integer value specifying the number of tubes across 
%         the width of the river. Default value = 20.
%       NoVertTubes = Integer value specifying the number of layers of 
%         tubes in each vertical. Default value = 5.
%       CellsPerXs = Integer value specifying the number of delft3D cells
%         between each streamtubes cross-section. Default value = 2.
%       TotalFlow = Optional user specifed total flow. If not supplied
%         calculated mean cross-section flow is used instead.
%       MaxDryCellsInTube = Optional user specified threshold number of dry 
%         cells across the width of any given streamtube. Streamtubes with 
%         more dry cells than this (i.e. streamtubes which span islands)
%         will be identified in the mask layer. Default value = 3.
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
%       Mask = cell array with size = [No_of_cross-sections,1] where each
%         cell is a matrix with size = [No_of_nodes,3] (same as Nodes).
%           Contains 1 or 0 values indicating active and inactive Nodes. 
%   Richard Measures, Gu Stecca, NIWA
%
%   See also streamtubeXS delft3d_multiflows2streamtubes

%% ############### TO DO ################
% - Improve calculation of depth at cell faces taking into account dpsopt
%   and dpuopt... confusing!

%% Get File names if not supplied
if (~exist('MdfFName','var')||isempty(MdfFName))
    [MdfFName,FilePath] = uigetfile('*.mdf','Select the delft3d model definition file');
    if isequal(MdfFName,0)
        error('User selected Cancel')
    end
    MdfFName = fullfile(FilePath,MdfFName);
    clear FilePath
end

if (~exist('StreamtubesFName','var')||isempty(StreamtubesFName))
    [StreamtubesFName,FilePath] = uiputfile('*.txt','Select file for streamtubes output', 'Streamtubes.txt');
    if isequal(StreamtubesFName,0)
        error('User selected Cancel')
    end
    StreamtubesFName = fullfile(FilePath,StreamtubesFName);
    clear FilePath
end

if (~exist('MaskFName','var')||isempty(MaskFName))
    [MaskFName,FilePath] = uiputfile('*.txt','Select file for streamtubes mask output', 'Mask.txt');
    if isequal(MaskFName,0)
        error('User selected Cancel')
    end
    MaskFName = fullfile(FilePath,MaskFName);
    clear FilePath
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

%% Set defaults if inputs not provided
if (~exist('NoHorizTubes','var')||isempty(NoHorizTubes))
    NoHorizTubes = 20;
end
if (~exist('NoVertTubes','var')||isempty(NoVertTubes))
    NoVertTubes = 5;
end
if (~exist('CellsPerXs','var')||isempty(CellsPerXs))
    CellsPerXs = 2;
end
if (~exist('ks','var')||isempty(ks))
    if isfield(MDF.keywords, 'ccofu') && ~isempty(MDF.keywords.ccofu) && isnumeric(MDF.keywords.ccofu)
            ks = (MDF.keywords.ccofu * 8.1 * 9.81^0.5)^6 * ones(MDF.keywords.mnkmax(2)-2, MDF.keywords.mnkmax(1)-2);
            disp('roughness value is read from .mdf file.')
    elseif isfield(MDF.keywords, 'filrgh')
            disp('.rgh file exists, roughness values are read from .rgh file.');

            % Split the ModelPath into its components
            modelPathComponents = split(ModelPath, filesep);
            % Split the relative path to .rgh file into its components
            relativePathComponents = split(MDF.keywords.filrgh, filesep);
            
            % Traverse the relative path and adjust the model path components
            for i = 1:length(relativePathComponents)
                if strcmp(relativePathComponents{i}, '..')
                    % Go up one directory level
                    modelPathComponents(end) = [];
                else
                    % Append the relative path component to the model path
                    modelPathComponents{end + 1} = relativePathComponents{i};
                end
            end
            % Construct the full path
            RghPath = fullfile(modelPathComponents{:});

            % Read the spatially varying roughness file
            Manning_rgf = dlmread(RghPath);
            if ismatrix(Manning_rgf) && size(Manning_rgf, 1) > 1 && size(Manning_rgf, 2) > 1
                disp('ks is a 2D matrix.');
            end
            a=1; b=1;
            for i = 1:(size(Manning_rgf, 1)/2) % we use only the u component of Manning
                for j = 1:size(Manning_rgf, 2)
                    if Manning_rgf(i,j) ==0
                        continue
                    end
                    ks(a,b) = Manning_rgf(i,j);
                    if mod(b,MDF.keywords.mnkmax(1)) == 0
                        a=a+1; b = 1;
                    else b=b+1;
                    end
                end
            end
            ks = (ks .* 8.1 * 9.81^0.5).^6;
            ks = ks(2:end-1, 2:end-1);
    else ks = 0.1 * ones(MDF.keywords.mnkmax(2)-2, MDF.keywords.mnkmax(1)-1);
        disp('setting ks to default value of 0.1m.');
    end
end
if (~exist('MaxDryCellsInTube','var')||isempty(MaxDryCellsInTube))
    MaxDryCellsInTube=3;
end

%% Read model results

% Set default output time if not provided
if (~exist('OutputTimeID','var')||isempty(OutputTimeID))
    OutputTimeID = T.nt_loaded;
end

% Get depth at specified time
S0 = vs_get(H,'map-const','DPS0',{2:Grid.nmax-1,2:Grid.mmax-1},'quiet');
WL = vs_get(H,'map-series',{OutputTimeID},'S1',{2:Grid.nmax-1,2:Grid.mmax-1},'quiet');

Depth = S0 + WL;

Depth(isnan(Grid.cen.y)) = NaN;
DepthThreshold = MDF.keywords.dryflc;
Depth(Depth<DepthThreshold) = 0;
clear S0 WL

% Remove disconnected ponds of water
WetConnected = bwareaopen(Depth>=DepthThreshold, 10000, 4);
Depth(WetConnected==0) = 0;
Depth(Depth<=0) = 0;

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
    ks = flipud(ks);
    Xcor  = flipud(Grid.cor.x);
    Ycor  = flipud(Grid.cor.y);
else
    % N direction dominant (so transpose, flipud not necessary)
    Vel   = V1';
    Depth = Depth';
    ks = ks';
    Xcor  = Grid.cor.x';
    Ycor  = Grid.cor.y';
end

if mean(mean(Vel)) < 0
    % flow in negative direction so fliplr to process from upstream to downstream
    Vel   = -fliplr(Vel);
    Depth = fliplr(Depth);
    ks = fliplr(ks);
    Xcor  = fliplr(Xcor);
    Ycor  = fliplr(Ycor);
end

% convert depths and roughness to cell faces
Depth = [Depth(:,1),Depth,Depth(:,end)];
Depth = (Depth(:,1:end-1) + Depth(:,2:end)) /2;
ks = [ks(:,1),ks,ks(:,end)];
ks = (ks(:,1:end-1) + ks(:,2:end)) /2;

% Tidy up
clear U1 V1

%% Build streamtubes

SelectedXs = [1:CellsPerXs:size(Depth,2)-1,size(Depth,2)]';
NoOfXs = size(SelectedXs,1);

Nodes = cell(NoOfXs,1);
Tubes = cell(NoOfXs,1);
Mask = cell(NoOfXs,1);
Flows = nan(NoOfXs,1);
% ks_Node= nan(NoOfXs,1);

% Loop through cross-sections
XsCount = 0;
for XsNo = SelectedXs'
    XsCount = XsCount+1;
    XS_X     = Xcor(:,XsNo);
    XS_Y     = Ycor(:,XsNo);
    XS_Vel   = Vel(:,XsNo);
    XS_Depth = Depth(:,XsNo);
    XS_ks = ks(:,XsNo);
  
    [Nodes{XsCount,1},Tubes{XsCount,1},Flows(XsCount,1), Mask{XsCount,1}, ks_Node{XsCount,1}] = ...
        streamtubeXS(XS_X, XS_Y, XS_Vel, XS_Depth, XS_ks, ...
                     NoHorizTubes, NoVertTubes, MaxDryCellsInTube);
end
clear XS_X XS_Y XS_Vel XS_Depth XS_Count Mask_Depth XS_Mask XS_ks

%% Write out streamtubes file

MeanFlow = mean(Flows);
if (~exist('TotalFlow','var')||isempty(TotalFlow))
    TotalFlow = MeanFlow;
else
    if (TotalFlow > (MeanFlow * 1.1)) || (TotalFlow < (MeanFlow * 0.9))
        warning('TotalFlow specified (%.3f m^3/s) is significantly different to meanflow at cross sections (%.3f m^3/s)',TotalFlow,MeanFlow)
    end
end

plotStreamtubes2d(Nodes,Tubes,Mask);

writeStreamtubes(Nodes, Tubes, StreamtubesFName, TotalFlow, ks_Node, Mask, MaskFName)

clear ks StreamtubesFName TotalFlow Flows MeanFlow ks_NMode MaskFName MDF Grid Depth Vel S0 WL

end

