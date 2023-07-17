# Lactic Acid Inhibition
## Simulation Scripts and Data
In the directory ```simulations/``` there is code for the optimization of a D-lactic acid production process with acid inhibition: ```mcp_v2**.jl```. In the respective directory ```mcp_v2**/``` contains the solutions of the simulations:
* ```summary.txt```, a text summary of the results
* ```df.csv```, the integrated process variables including their rates and fluxes.
* ```q.csv```, the fluxes of the FBA solution at every FE.
* ```vlb.csv```and ```vub.csv``` the flux bounds of the FBA model.
* ```variables.jld2``` a collection of all process parameters in hdf5/jld2 format.
   
Table of Simulations:

|ID  | Description|
|----|----|
mcp_v200 | Optimize final product mass, multi stage
mcp_v205  | Optimize final product mass over final time, multi-stage
mcp_v232  | Optimize final product mass over final time, two-stage
mcp_v233  | Optimize final product mass over final time, single-stage

## Analysis Scripts
The Notebook ```Plot_Simulations.ipynb``` plots the process variables and the production envelope of the process.
There is also a comparison of one-, two- and multi-stage processes.

## FBA Model Files and Scripts
```iJO1366.xml```, ```lbi.csv```, ```ubi.csv```, ```Si.csv```, and ```postprocessing_02.jl``` contains the FBA model (bounds) and some helper scripts.


