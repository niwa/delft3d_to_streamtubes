%% Run delft3D streamtubes function

%% Add Open Earth Matlab Toolbox 
% addpath C:\Projects\OpenEarth\matlab\
% oetsettings ('quiet')

%% Reading data
%Mdf file names
MdfNames = 'C:\Users\torokg\OneDrive - NIWA\projects\ENS24505_Waiau\model\discharges\mdf_files.dat';
    % Open the file for reading
    fileID = fopen(MdfNames, 'r');
    % Check if the file opened successfully
    if fileID == -1
        error('Failed to open the file.');
    end
    fileData = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    fileNames = fileData{1};

% Reading discharge values
discharges = dlmread('C:\Users\torokg\OneDrive - NIWA\projects\ENS24505_Waiau\model\discharges\discharges.dat');

%% Specify settings and Generate streamtubes
for FlowNumber = 1:length(fileNames)
%     FlowNumber = 1;
    MdfFName = fileNames{FlowNumber};
    StreamtubesFName = sprintf('Streamtube_constks_Q%g.txt', discharges(FlowNumber));
    MaskFName = sprintf('Mask_constks_Q%g.txt', discharges(FlowNumber));
    NoHorizTubes = 20;
    NoVertTubes = 5;
    CellsPerXs = 2;
    TotalFlow = discharges(FlowNumber);
    OutputTimeID = [];
    ks = [];
    MaxDryCellsInTube = 3;

    [Nodes, Tubes, Mask] = delft3d_streamtubes(MdfFName, OutputTimeID, StreamtubesFName, MaskFName, ks, NoHorizTubes, NoVertTubes, CellsPerXs, TotalFlow);
end

    %% Plot test output for a single flow
%    Be carefull - if there are lots of tubes/cross-sections a 3d plot can 
%    run out of memory.

%     3D plot:
%     plotStreamtubes(Nodes,Tubes);
%     2D plot:
%     plotStreamtubes2d(Nodes,Tubes,Mask);
% end

