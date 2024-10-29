# Plos One Submission - Core Code Repo

The code contained within this repo provides function for 
(1) ODE simulation of fluorohpore population under single pulse excitation
(2) Stochastic simulation of fluorophore populations by employing Gillespie Algorithm (GA) under single pulse excitation

The core function for the ODE simulation is named ode_sol.m 

The core function for the GA simulation is named ssa_engine.m 

Each function takes in arguements about the simulation duration, state rates, fluorophore photophysical properties, and excitation properties. These core codes were used to generate the time and state vectors decsribing the temporal evolution of the fluorohpore populations. With these vector infomration, one can analyze a variety of metrics such as average rise time to steady state, steady-state emission photon flux, quantum state occupany levels, and generating the power spectral density of the photon arrival (by apporximating the PSD from the photon arrival times). 
