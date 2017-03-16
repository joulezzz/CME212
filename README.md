# CME212 Submission Repository

This is your CME212 homework repository. Feel free to edit your README and any other files!

I added 4 cores in the setting in virtual box to allow the Virtual Linux System on my mac to run the mass_spring executable on potentially 4 cores. The following results are shown below in a table after running on different grid sizes for grid0, grid1, grid2, and grid3. 

Table of Results
 
|                                                | thrust::system::omp::par | thrust::system::detail::sequential::seq |
|------------------------------------------------|--------------------------|-----------------------------------------|
| ./mass_spring data/grid0.nodes data/grid0.tets | 148.637                  | 155.018                                 |
| ./mass_spring data/grid1.nodes data/grid1.tets | 141.23                   | 178.846                                 |
| ./mass_spring data/grid2.nodes data/grid2.tets | 292.818                  | 139.668                                 |
| ./mass_spring data/grid3.nodes data/grid3.tets | 419.49                   | 322.858                                 |