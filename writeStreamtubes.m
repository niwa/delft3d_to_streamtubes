function writeStreamtubes(Nodes, Tubes, OutFName)
%WRITESTREAMTUBES output streamtubes to text file
%   writeStreamtubes(Nodes, Tubes, OutFName)
%
%   Richard Measures 2016
%
%   See also delft3d_streamtubes

%% Open file for write access
FID=fopen(OutFName,'w');

%% Write header info
fprintf(FID, 'number of x-secs: %i\n', size(Nodes,2));
fprintf(FID, 'number of flow tubes: %i\n', );


end

