%% Run delft3D streamtubes function

%% Add Open Earth Matlab Toolbox 
addpath C:\Projects\OpenEarth\matlab\
oetsettings ('quiet')

%% Specify settings
MdfFName = 'C:\MyDirectory\MyDelft3DModelDefinitionFile.mdf';
StreamtubesFName = 'MyNewStreamtubes.txt';
ks = 0.1;
NoHorizTubes = 20;
NoVertTubes = 5;
CellsPerXs = 2;

%% Generate streamtubes
[Nodes, Tubes] = delft3d_multiflows2streamtubes(MdfFName, StreamtubesFName, ks, NoHorizTubes, NoVertTubes, CellsPerXs);

%% Plot test output for a single flow
%    Be carefull - if there are lots of tubes/cross-sections a 3d plot can 
%    run out of memory.
FlowNo = 11;

% 3D plot:
plotStreamtubes(Nodes{FlowNo,1},Tubes{FlowNo,1});
% 2D plot:
plotStreamtubes2d(Nodes{FlowNo,1},Tubes{FlowNo,1});
