% Lukas WinklerPrins
% lukas_wp@berkeley.edu

% Last edited 30 April 2020

close all
clear all

for ss = 1:4
    sensor_choice = ss;
    make_contour_plot
end

[~,si] = min(abs(start_time - datenum(buoy_spectra.time)));
[~,ei] = min(abs(end_time - datenum(buoy_spectra.time)));
mean_buoy_spectra = mean(buoy_spectra.spectra(si:ei,:));

Cg_buoy = (g/(4*pi)).*(1./buoy_spectra.freq); % Using deep-water assumption

buoy_energy_flux = mean_buoy_spectra.*Cg_buoy;

figure
loglog(freq,energy_flux_1);
hold on
loglog(freq,energy_flux_2);
loglog(freq,energy_flux_3);
loglog(freq,energy_flux_4);
loglog(buoy_spectra.freq,buoy_energy_flux);
xlabel('frequency (Hz)');
ylabel('Energy Flux (m^3/s)');
title(['Energy Flux from ' datestr(start_time) ' to ' datestr(end_time) ' (high local winds)']);
legend('LL','SB','PN','TP','NDBC Buoy') % For 4th Deployment

dEF_LL_to_SB = energy_flux_2 - energy_flux_1;
dEF_SB_to_PN = energy_flux_3 - energy_flux_2;
dEF_PN_to_TP = energy_flux_4 - energy_flux_3;

figure
semilogx(freq,dEF_LL_to_SB);
hold on
semilogx(freq,dEF_SB_to_PN);
semilogx(freq,dEF_PN_to_TP);
xlabel('frequency (Hz)');
ylabel('Energy Flux Difference (m^3/s)');
title(['Change in Energy Flux Going Into the Bay, ' datestr(start_time) ' to ' datestr(end_time)]);
legend('LL to SB','SB to PN','PN to TP');
