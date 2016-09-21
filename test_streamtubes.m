%% Test delft3D streamtubes function

%% Add Open Earth Matlab Toolbox 
addpath C:\Projects\OpenEarth\matlab\
oetsettings ('quiet')

%% Specify settings
%MdfFName = 'C:\Projects\Research\SWAP_3Rivers\UpperCust\UpCust_refine.mdf';
MdfFName = 'S:\measuresrj\Oreti\Oreti_RampedQ_highflow.mdf';
StreamtubesFName = 'streamtubes.txt';
ks = 0.1;
NoHorizTubes = 10;
NoVertTubes = 5;
CellsPerXs = 20;

%% Generate streamtubes
[Nodes, Tubes] = delft3d_multiflows2streamtubes(MdfFName, StreamtubesFName, ks, NoHorizTubes, NoVertTubes, CellsPerXs);

%% Plot test output for a single flow
FlowNo = 2;
plotStreamtubes(Nodes{FlowNo,1},Tubes{FlowNo,1});
