% Lukas WinklerPrins
% lukas_wp@berkeley.edu
% Last Edited 16 July 2020

% First, need to find energy flux for each frequency bucket. 
% Energy Flux = Group Celerity * Energy Density 
% (really needs to then be multiplied by rho*g*(1/8)...

% H = 0.5; % m

% E = 0.5*rho*g*(H/2).^2;

% u_b = ((H/2).*W)./sinh(kh);
% 
% figure
% plot(freq,Cp)
% hold on
% plot(freq,Cg)
% xlabel('Frequency (Hz)');
% ylabel('Speed (m/s)');
% legend('C_p','C_g');

% epsilon_f = rho*c_f*((W./sinh(kh)).^3).*((H.^3)./(6*pi));
% epsilon_f = rho*c_f*u_b.^3;
% linear_epsilon_f = epsilon_f./Cg;