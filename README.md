# LocalizedFlow_Pockmarks
Model code developed in "Hydromechanical Simulations Reveal the Climatic-Relevance of CO2 Release from Pockmark on the Chatham Rise-Bounty Trough, New Zealand"
Katrina C. Magno, J. I. T. Hillman, L. Räss, L. Stott, J. Suckale (in prep. 2025)
This collection of code includes a (1) 2D hydromechanical model (Julia) that simulates the formation of localized fluid pathways via nonlinear percolation of fluid through a viscous, deforming sediment matrix and (2) simulation data analyses code (MATLAB) that identifies pockmarks as localized fluid pathways that cross an idealized seafloor and outputs vertical flux [km/yr] vs area [km^2] plots and estimates a range of fluid volume transport rates by scaling according to observational pockmark size data from Chatham Rise, New Zealand.

## 2D Hydromechanical Model (Julia):
### HydroMech2D_Main.jl:
Main code to run HM model
  - Saves out vertical flux (qDy), time, and porosity (ϕ)
  - Adjust HM Model GPU capability via "const USE_GPU". Setting to "true" enables HydroMech2D_Main.jl to be run on GPUs. Model was submitted using sbatch to the Stanford Sherlock HPC       Cluster

### HydroMechDataTypes.jl:
  - Specify poromechanical parameters and initial conditions
  - For distributed flow: R = 1.0 and n_k = 1.0 (Figure 2a-d in Magno et al., 2025)
  - For localized flow: R > 1.0 and n_k > 1.0 (Figure 2e-h in Magno et al., 2025)

### HydroMechFunctions.jl:
  - Functions required to solve coupled Darcy + Stokes flow equations

## Simulation Data Analyses (MATLAB):
### Gen_FluxvsArea.m:
Inputs: vertical flux data, time, and porosity CSV files from HM simulation. (Optional: Simulation data included in ... to reproduce Figure 3 d-e)
Outputs: Plots of Average Vertical Flux/pockmark [km/yr] vs pockmark area [km^2] and range of the fluid volume transport rates [km^3]/yr from 476 Chatham Rise pockmarks. Assumes        fluid is 100% CO2 and is in either liquid or gaseous phase
Initial Flow:
  - Use "natsort.m" to sort input CSV files according to simulation time (Stephen23 (2025). Customizable Natural-Order Sort (https://www.mathworks.com/matlabcentral/fileexchange/34464- 
    customizable-natural-order-sort), MATLAB Central File Exchange.)
  - For each n_k value, open qDy (M(vertical model resolution x horizontal model resolution)), Phi (M(vertical model resolution x horizontal model resolution)), and time (vector)           files. The number of files corresponds to the number of simulation time steps
  - For each time step/file, collect data from row (horz_slice) in qDy and Phi that represents the idealized seafloor 200 m above the fluid reservoir, creating new matrices qdy_slices      M(#time steps x horizontal model resolution) and phi_slices M(# time steps x horizontal model resolution). These matrices describes the temporal evolution of vertical fluid flux        and porosity along the idealized seafloor
  - In order to supress background noise, threshold data by setting flux values lower than 2 (the background) to zero
    
Identifying Pockmarks. For each time step/row in qdy_slices and phi_slices:
  - Identify peaks in row of qdy_slices by plotting the vertical flux data against horiztonal location. Peaks in vertical flux translate to the formation of pockmarks at the                idealized seafloor. Use MATLAB Signal Processing Toolbox function "findpeaks" to identify peak locations, prominences, and full width at half maximums (FWHM)
  - Estimates the pockmark area by using the FWHM/2 as the radius of a circular pockmark
  - Using location and FWHM, isolate the entire curve of each peak and spatially average porosity and vertical flux for each pockmark at that time step
  - Fit each curve to a Gaussian distribution. Using MATLAB "fit" and "gauss1"
  - Keep track of locations of peaks at each time step and compare location to previously recorded locations. More specficially, if new location x_oi is equal to or +/- 1 away from a       previously recorded location x_01, assign the vertical flux, porosity, and pockmark area estimates to x_01. This is done to achieve a degree of uniqueness in our final pockmark         estimates
  - Using "calc_vol.m" included in this repository, expand each Gaussian distribution to a 2D Gaussian surface
  - Estimate the fluid volume transport rate [km^3/yr] by multiplying the integral [km^2] of the 2D Gaussian by the peak vertical flux value [km/yr] for each pockmark
  - Scale volume transport rate by mean porosity
    
After each row of qdy_slices is processed, time-average the spatially-averaged vertical flux and pockmark areas to get the average vertical flux per pockmark [km/yr] and associates average pockmark area [km^2]. Generate plots as Figure 3 d-e. Also time-average the estimate fluid volume transport rates to get average fluid volume transport rates associated with average pockmark area.

Estimate fluid volume transport rates for 476 Pockmarks on Chatham Rise, New Zealand:
  - Select the simulated pockmarks that correspond to the same size classes of the 476 pockmarks observed on the Chatham Rise
  - Scale simulated fluid volume transport rates by the size distribution of 476 pockmarks and sum quantities to get the total fluid volume transport rate [km^3/yr] from the         
    simultaneous eruption of fluid from the 476 pockmarks
  - Assume fluid is 100% CO2 and is in either liquid or gaseous phase to get mass transport rate estimates of CO2 [PgCO2/yr]

### DOI: https://zenodo.org/doi/10.5281/zenodo.13328488
