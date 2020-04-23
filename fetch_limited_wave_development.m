function [Hm0,Tm] = fetch_limited_wave_development(F,U)
% FETCH_LIMITED_WAVE_DEVELOPMENT returns H_m0 (significant wave height) and
% Tm (significant period) based on Shore Protection Manual's equations 3-33
% and 3-34, for fetch-limited (rather than duration-limited) conditions. F
% is a fetch length (in meters) and U is the wind speed (in m/s), scalars.

% Lukas WinklerPrins
% lukas_wp@berkeley.edu
% Last Edited 22 April 2020

U_A = 0.71*(U)^1.23;

g = 9.80665; % m/s2

Hm0 = (U_A^2)*(1/g)*1.6*(10^-3)*(g*F/(U_A^2))^(1/2);
Tm =  (U_A/g)*2.857*(10^-1)*(g*F/(U_A^2))^(1/3);

end