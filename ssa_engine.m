%%This code was modified from the original source code written by Nezar Abdennur, 2012 

%%Author: Nick Vitale 

%Function inputs: 
% 1. kl - 1/tau of fluorophore lifetime 
% 2. rx - rate of fluorophore promotion from ground to excited state 
% 3. QE - qantum efficiency of fluorophore
% 4. duration - max duration of reaction to observe under CW conditions.
% NOTE: At the end of the specified duration, it is assumed the excitation
% source deactivates with ideal fall time. 
% 5. end_time - end time of the simulation. Ensure that this number exceeds
% the time one wants to observe the TG portion of the waveform. 
% 6. nf - number of fluorophores to simulate at once 

%Outputs
%t_photon - times of photon emission 
%t - time steps of the GA simulation 
%x - resulting state matrix of the GA simulation 

function [t_photon,t,x] = ssa_engine(kl,rx,QE,duration,end_time,nf)

%% Rate constants
p.relax = kl; %Rate at which fluorophores relax from excited to ground state (includes both radiative and non-radiative processes)
p.excited = rx; %Rate at which fluorohpores are promoted from the ground to excited state
p.t_dur = duration;
p.t_end = end_time;
p.QE = QE;

%% Initial state
tspan = [0,p.t_dur+p.t_end,p.t_dur];
x0    = [nf, 0];  %Initialize the stoich. matrix - [0]: # of fluorophores in ground state, [1]: # of fluorohpores in excited state 

%% Specify reaction network
pfun = @propensities_2state;
stoich_matrix = [ -1  1    %fluorophore removed from ground and sent to excited
                   1  -1 ]; %fluorophore removed from excited to ground

%% Run simulation
%Last entry is set as the max number of reactions allowed to occur.
%Increase if needed. 
[t_photon,t,x] = directMethod_mod(stoich_matrix, pfun, tspan, x0, p,[],100E6);

end

%Propensity network for ground and excited states. You can expand this
%linearly given that the states and their rates are known. 

function a = propensities_2state(x, p)

ng   = x(1); %Number of fluorophores in the ground state 
ne   = x(2); %Number of fluorophores in the excited state 

a = [p.excited*ng;            %ground to excited
     p.relax*ne;];            %excited to ground
end
