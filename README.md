# LocalizedFlow_Pockmarks
This collection of code includes a (1) 2D hydromechanical model (Julia) that simulates the formation of localized fluid pathways via nonlinear percolation of fluid through a viscous, deforming sediment matrix and (2) post-simulation code (MATLAB) that identifies pockmarks as localized fluid pathways that cross an idealized seafloor and outputs vertical flux [km/yr] vs area [km^2] plots and estimates a range of fluid volume transport rates by scaling according to observational pockmark size data from the Chatham Rise, New Zealand.

# 2D Hydromechanical Model (Julia):
HydroMech2D_Main.jl:
  - Main code to run HM model and save out vertical flux (qDy), time, and reference porosity (ϕ0) 
HydroMechDataTypes.jl:
  - Specify poromechanical parameters and initial conditions
HydroMechFunctions.jl:
  - Functions required to solve coupled Darcy + Stokes flow equations

# Post-Simulation Analysis (MATLAB):
Gen_FluxvsArea.m:
  - Identifies pockmarks in horizontal slices of flux data by fitting Gaussian curves to peaks in the slices
  - Estimates the pockmark area by using the FWHM as the radius of a circular pockmark
  - Estimates the average porosity of a given pockmark size
  - Estimates the average vertical flux [km/yr] from a pockmark of a given size
  - Links simulations to observational data by scaling simulated quantities by the size distribution of 476 pockmarks observed on the Chatham Rise, NZ
  - Finally, generates a range of the fluid volume transport rate [km^3]/yr from 476 Chatham Rise pockmarks. Assumes fluid is 30-70% CO2 and is in either liquid or gaseous phase
