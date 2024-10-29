%This function was directly modified from the GA algo originally written by Nezar Abdennur, 2012


function [t_photon, t, x] = directMethod_mod( stoich_matrix, propensity_fcn, tspan, x0,...
                                  params, output_fcn, MAX_OUTPUT_LENGTH)
%DIRECTMETHOD Implementation of the Direct Method variant of the Gillespie algorithm
%   Usage:
%       [t, x] = directMethod( stoich_matrix, propensity_fcn, tspan, x0 )
%       [t, x] = directMethod( stoich_matrix, propensity_fcn, tspan, x0, params )
%       [t, x] = directMethod( stoich_matrix, propensity_fcn, tspan, x0, params, output_fcn )
%
%   Returns:
%       t:              time vector          (Nreaction_events x 1)
%       x:              species amounts      (Nreaction_events x Nspecies)    
%
%   Required:
%       tspan:          Initial and final times, [t_init, t_final].
%
%       x0:             Initial species amounts, [S1_0, S2_0, ... ].
%
%       stoich_matrix:  Matrix of stoichiometries (Nreactions x Nspecies).
%                       Each row gives the stoichiometry of a reaction.
%
%       prop_fcn:       Function that calculates reaction propensities.
%                       Function should be of the form
%                           ac = f( xc, params )
%                       where xc is the current state [S1, S2, ...], and
%                       params is the user-defined rate parameters.
%                       The function should return column vector ac of 
%                       propensities (Nreactions x 1) in the same order as
%                       the reactions given in stoich_matrix.
%
%   Optional:
%       params:         User-defined parameters, passed to prop_fun (e.g.,
%                       a struct of rate constants) <default=[]>
%
%       output_fcn:     Arbitrary function with signature
%                           status = f( tc, xc )
%                       The output_fcn is passed the current time and state
%                       after each step of the simulation. It can be used
%                       to locate events, monitor progress, write data,
%                       etc. If it returns 1, the simulation terminates.
%                       <default=none>
%
%   Reference: 
%       Gillespie, D.T. (1977) Exact Stochastic Simulation of Coupled
%       Chemical Reactions. J Phys Chem, 81:25, 2340-2361.
%


if ~exist('MAX_OUTPUT_LENGTH','var')
    MAX_OUTPUT_LENGTH = 1000000;
end
if ~exist('output_fcn', 'var')
    output_fcn = [];
end
if ~exist('params', 'var')
    params = [];
end
    
%% Initialize
num_species = size(stoich_matrix, 2);
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, num_species);

%%Pre-initialize random numbers ahead of time
r1 = rand(1,MAX_OUTPUT_LENGTH);
r2 = rand(1,MAX_OUTPUT_LENGTH);
r3 = binornd(1,params.QE,1,MAX_OUTPUT_LENGTH);

T(1)     = tspan(1);
X(1,:)   = x0;

%Initialize arrays to count photon arrival times 
rxn_count = 1;
t_photon = [];
photons = 0;
photon_count=[];
FLAG = 1;

%% MAIN LOOP
while T(rxn_count) < tspan(2)        
    
    % Calculate reaction propensities
    a = propensity_fcn(X(rxn_count,:), params);
    
    % Sample earliest time-to-fire (tau)
    a0 = sum(a);
    tau = -log(r1(rxn_count))/a0; 
    
    % Sample identity of earliest reaction type to occur (mu: 1= ground to excited, 2 = excited to ground)
    mu = find((cumsum(a) >= r2(rxn_count)*a0), 1,'first');

    %If the number of reactions that have occured exceeds the allocated
    %number, terminate
    if rxn_count + 1 > MAX_OUTPUT_LENGTH
        t = T(1:rxn_count);
        x = X(1:rxn_count,:);
        warning('SSA:ExceededCapacity',...
                'Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
        return;
    end

    %If the reaction type is excited to ground, then we choose a random
    %y/n decision based on a binomial distribution determined by the QE of
    %the fluorohpore. If photon fired, store the time.
    if(mu == 2)
        if(r3(rxn_count) == 1)
            t_photon= [t_photon,T(rxn_count)];
            photons = photons+1;
        end
    end


    %If the next reaction occurs at the time past the duration of the pulse
    %steady state, then we need to stop, update the rate of ground to
    %excited to be 0, and then resume simulation to catch the decay of the
    %probe population
    if((FLAG)&&(T(rxn_count)   + tau) > tspan(3))
        T(rxn_count+1) = tspan(3);
        X(rxn_count+1,:) = X(rxn_count,:) + [0,0];
        photon_count(rxn_count+1) = photons;
        params.excited = 0;
        FLAG=0;
    else
    %Update and keep moving
        T(rxn_count+1)   = T(rxn_count)   + tau;
        X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);
        photon_count(rxn_count+1) = photons;
    end

    fprintf('Sim Time is now: %d\n',T(rxn_count));

    rxn_count = rxn_count + 1;
    


    if ~isempty(output_fcn)
        stop_signal = feval(output_fcn, T(rxn_count), X(rxn_count,:)');
        if stop_signal
            t = T(1:rxn_count);
            x = X(1:rxn_count,:);
            warning('SSA:TerminalEvent',...
                    'Simulation was terminated by OutputFcn.');
            return;
        end 
    end


   
end  

% Return simulation time information of the photon arrivals
t = T(1:rxn_count);
x = [X(1:rxn_count,:),photon_count(1:rxn_count)'];

if t(end) > tspan(2)
    t(end) = tspan(2);
    x(end,:) = [X(rxn_count-1,:),photon_count(rxn_count-1)'];
end    

end

