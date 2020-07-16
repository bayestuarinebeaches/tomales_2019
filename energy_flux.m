% Lukas WinklerPrins
% lukas_wp@berkeley.edu
% Last Edited 16 July 2020

% First, need to find energy flux for each frequency bucket. 
% Energy Flux = Group Celerity * Energy Density 
% (really needs to then be multiplied by rho*g*(1/8)...

g = 9.80665;

% All of these are vectors, values for each frequency (really, W) value
h = mean(trimmed_depth_signal);
% h = 4;

kh = qkhfs(W,h); % Oy, taking in a lot of time here...
k = kh./h; % a bit sloppy
Cp = ((g./k).*tanh(kh)).^0.5;
Cg = Cp.*(0.5+kh./sinh(2.*kh));

c_f = 0.01; % Per Tucker & Pitt p321
rho = 1018; % arbitrary guess
H = 0.5; % m

E = 0.5*rho*g*(H/2).^2;

u_b = ((H/2).*W)./sinh(kh);

% epsilon_f = rho*c_f*((W./sinh(kh)).^3).*((H.^3)./(6*pi));
epsilon_f = rho*c_f*u_b.^3;
% linear_epsilon_f = epsilon_f./Cg;