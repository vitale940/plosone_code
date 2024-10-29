%Function used to simulate via ODE solver
% Inputs: 
% nf - number of fluorophores to be simulated 
% ss_dur - total time to observed rise + steady state regime of waveform 
% sim_stop - stop time of the simulation. This parameter should > ss_dur in
% order to observe the decay of the emission after pulse deactivation 
%illum_surface - surface (in cm^2) of the illuminated surface to compute
%excitation flux 
%probe_lifetime - lifetime of fluorophore in seconds 
%bleach_lifetime - lifetime of bleaching process in seconds 
%quant_eff - quantum efficiency of fluorohpore 
%molar_abs - molar absorption coefficient (in mol^-1 cm^-1)
%exc_power - power intensity of the excitation source in watts 
%excitation_wavelength - center wavelength of excitation source 
%probe_surface - surface on which probes are spotted. Generally this is set
%equal to illum_surface


%Outputs 
%t_ode - Output time vector from ODE solution
%xa_ode - Output vector containting the temporal evolution of the system.
%       - [1] = ground state solution
%       - [2] = excited state solution
%       - [3] = bleached state solution
%       - [4] = emission flux solution

function [t_ode, xa_ode] = ode_sol(nf, ss_dur, sim_stop, illum_surface,probe_lifetime, bleach_lifetime, quant_eff, molar_abs, exc_power, excitation_wavelength, probe_surface )

    %% Excitation Source Timing, Power, Surface Area Information, and Signal Constrcution
    %Initial Conditions
    %nf is initilized via input
    ne = 0;
    nb = 0;
    nph = 0;
    
    %Timing
    start = 0;
    dur = ss_dur; %Duration of the pulse (seconds)
    rise = 1E-12; %RC Rise Time of pulse (seconds)
    fall = 1E-12; %RC Fall Time of pulse (seconds)
    stop = sim_stop;
    delay = 0;
    
    %Illumination Surface Area
    S1 = illum_surface; %in cm^2
    
    %% Fluorophore Information
    taul =probe_lifetime; % Probe lifetime (seconds)
    kl = 1/taul;
    taub = bleach_lifetime; % Bleaching lifetime (seconds)
    kb = 1/taub;
    QE = quant_eff; % Quantum Efficiency
    epsilon = molar_abs; %Molar Absorption Coefficient (1/M*cm)
    
    %Power and Wavelength
    power = exc_power; %Power of source (watts)
    cw_wavelength = excitation_wavelength; %Center wavelength of light (um)
    
    %% ODE Simulator (ODE 15s)
    clear t_ode xa_ode
    tic(); %Start timer
    
    %System Parameter Infromation
    %x(1) = ng - Number of ground state fluorophores
    %x(2) = ne - Number of excited state fluorophores
    %x(3) = nb - Number of bleached state fluorophores
    %x(4) = nph - Number of photons generated
    
    S = probe_surface; %Surface area of interest on which the fluorophores are clustered.(cm^2) ##Provided as separate parameters in case surface area you are illuminating is different than surface area used to compute the source flux.
    
    opts = odeset('Refine',8); %Simulation accuracy setting. Forces ODE15s to get close to DOE45 accuracy.
    init = [nf ne nb nph]; %Initial paramters [ng,ne,nb,nph]
    
    %Computation of the flux based on the provided power and surface area 
    Eev = 1.24/cw_wavelength; %Energy of single photon at specified wavelength (eV)
    Ej = Eev*1.602E-19; %Energy of single photon at specified wavelength (joules)
    num_photons = power/Ej; %Rate o f photons from source (photons/s)
    flux = num_photons/S1;  %Computed Photonic Flux (photons/cm^2*s)   
    
    [t_ode,xa_ode] = ode15s(@(t,y) fun(t,y,flux,taul,taub,S,epsilon,QE,dur,delay),[start stop],init);
    
    r2 = toc(); %Stop timer
    
    fprintf('Simulation Time: %f Seconds\n',r2)

end
