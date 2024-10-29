%Function used for simulating ODE 
% Inputs: 
% t - time variable allocated for the ODE 
% y - variable that we used to represent the ODE entries 
% IT - vector of simulated LED source (generatored from
% light_source_constructor) 
% taul - lifetime of probe (seconds) 
% taub - bleaching lifetime of probe (seconds)
% surface - surface area where fluorophores are spotted (cm^2) 
% molar_co - molar absorpotion coefficient (1/M*cm) 
% eff - quantum efficieny of the fluorophore 

%Outputs 
%da - vector of the solutions wrt y 

function da = fun(t,y,IT,taul,taub,surface,molar_co,eff,dur,delay)

%Rate Constants 
kl = 1/taul; %Lifetime rate constant (1/s) 
kb =  1/taub; %Bleaching rate constant (1/s) 
QE = eff; %QE of probe 

%System information 
S = surface; %Surface area (cm^2)
M = 2.303*molar_co; %Adjusted Molar Coeff for the fact we are using the exponential form 
Av = 6.02E23; %Avagadros number 

%Using ODE45 to solve first (will get to close form later) 
%y(1) = ng - ground state probes 
%y(2) = ne - excited state probes 
%y(3) = nb - bleached probes
%y(4) = nph - photons generated 

%Interpolate the LED source vector 
%I = interp1(time,IT,t);
% 
if((t >=(delay)) && (t <= (dur+delay)))
    I = IT;
else
    I=0;
end
%ODE System 
da1 = kl.*y(2)-S.*I.*(1-exp(-M.*1000.*y(1)./(Av.*S))); %ng
da2 = S.*I.*(1-exp(-M.*1000.*y(1)./(Av.*S)))-kl.*y(2)-kb.*y(2); %ne
da3 = kb.*y(2); %nb
da4 = QE.*(kl+kb).*y(2);    %nph

%Stitch them into array to pass off to ODE solver 
da = [da1;da2;da3;da4];
end