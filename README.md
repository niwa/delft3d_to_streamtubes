# delft3d_to_streamtubes
Matlab functions for converting Delft3D model output to the streamtubes format 
required for input into the Cawthron NREI model. This software was developed 
with funding from NIWA's Environmental Flows research program.

##Contents:
The main streamtubes program:
* `delft3d_streamtubes.m`
  main function for generating streamtubes from delft3d results
* `streamtubeXS.m`
  function called in delft3d_streamtubes.m to generate streamtubes for individual 
  cross-sections
* `writeStreamtubes.m` 
  function called in delft3d_streamtubes.m to output streamtubes to text file
     
Wrapper to loop streamtubes generation for many flows:
* `delft3d_multiflows2streamtubes.m`
  identifies outputs from the end of stationary flow periods in an simulation with 
  a stepped inflow and applies delft3d_streamtubes.m to generate streamtube outputs 
  at each flow

Plotting finctions used to test streamtubes generation
 * `plotStreamtubes2d.m`
 * `plotStreamtubes.m`
 
Example applying functions to model results
 * `example_streamtubes.m`
 
See comments in individual files for more info.

Richard Measures, NIWA