%% Run delft3D streamtubes function

%% Add Open Earth Matlab Toolbox 
addpath C:\Projects\OpenEarth\matlab\
oetsettings ('quiet')

%% Specify settings
%MdfFName = 'F:\Delft3dProjects\Oreti\Oreti_RampedQ.mdf';
MdfFName = 'F:\Delft3dProjects\Oreti\Oreti_RampedQ_highflow.mdf';

StreamtubesFName = 'OretiStreamtubes.txt';
ks = 0.1;
NoHorizTubes = 20;
NoVertTubes = 5;
CellsPerXs = 2;

%% Generate streamtubes
[Nodes, Tubes] = delft3d_multiflows2streamtubes(MdfFName, StreamtubesFName, ks, NoHorizTubes, NoVertTubes, CellsPerXs);

%% Plot test output for a single flow
%    Be carefull - if there are lots of tubes/cross-sections a 3d plot can 
%    run out of memory.
FlowNo = 2;
%plotStreamtubes(Nodes{FlowNo,1},Tubes{FlowNo,1});
plotStreamtubes2d(Nodes{FlowNo,1},Tubes{FlowNo,1});
