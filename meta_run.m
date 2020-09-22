% Lukas WinklerPrins
% lukas_wp@berkeley.edu

% Last edited 17 July 2020

fprintf('Make sure you turned off sensor_choice, Lukas!!\n');

close all
clear all

run date_controls.m
run controls.m

n_sensors = length(cmab);

for ss = 1:n_sensors
    sensor_choice = ss;
    run_analysis
end

[~,si] = min(abs(start_time - datenum(buoy_spectra.time)));
[~,ei] = min(abs(end_time - datenum(buoy_spectra.time)));
mean_buoy_spectra = mean(buoy_spectra.spectra(si:ei,:));

Cg_buoy = (g/(4*pi)).*(1./buoy_spectra.freq); % Using deep-water assumption

buoy_energy_flux = mean_buoy_spectra.*Cg_buoy;

% figure
% loglog(freq,energy_flux_1);
% for ss = 2:n_sensors
%     hold on
%     eval(['loglog(freq,energy_flux_' num2str(ss) ');']);
% end
% loglog(buoy_spectra.freq,buoy_energy_flux);
% xlabel('frequency (Hz)');
% ylabel('Energy Flux (m^3/s)');
% title(['Energy Flux from ' datestr(start_time) ' to ' datestr(end_time) '']);
% legend(labels);
% % And last curve is NDBC Buoy

figure
semilogx(freq(n_leakage_ignore:end),running_S_1(n_leakage_ignore:end).*freq(n_leakage_ignore:end));
for ss = 2:n_sensors
    hold on
    eval(['semilogx(freq(n_leakage_ignore:end),running_S_' num2str(ss) '(n_leakage_ignore:end).*freq(n_leakage_ignore:end));']);
end
% semilogx(buoy_spectra.freq,mean_buoy_spectra.*buoy_spectra.freq);
xlabel('Frequency (Hz)');
ylabel('Energy Density (m^2 / Hz)');
title(['Energy Density from ' datestr(start_time) ' to ' datestr(end_time) ' using ' num2str(ea_spacing) '/' num2str(window_length) '/' num2str(instance_length) '.']);
legend(labels);
ylim([0 max_varpreserv_power]);
% And last curve is NDBC Buoy

figure
fg(1) = subplot(5,1,1);
for ss = 1:n_sensors
    hold on
    eval(['plot(window_times_' num2str(ss) ',Hs_wind_' num2str(ss) ');']);
end
xlabel('Datetime');
ylabel('H_s from Wind');
yyaxis right
plot(trimmed_times,trimmed_depth_signal,'k:');
legend(labels);

fg(2) = subplot(5,1,2);
for ss = 1:n_sensors
    hold on
    eval(['plot(window_times_' num2str(ss) ',Hs_swell_' num2str(ss) ');']);
end
xlabel('Datetime');
ylabel('H_s from Swell');
legend(labels);

fg(3) = subplot(5,1,3);
for ss = 1:n_sensors
    hold on
    eval(['plot(window_times_' num2str(ss) ',Hs_igw_' num2str(ss) ');']);
end
xlabel('Datetime');
ylabel('H_s from IGW');
legend(labels);

fg(4) = subplot(5,1,4);
yyaxis right
scatter(tbb_wind.time,tbb_wind.spd,'.');
%         scatter(datenum(wind.time),wind.spd,'.');
ylim([0 20]);
ylabel('Wind Speed (m/s)');
yyaxis left
hold on
scatter(tbb_wind.time,tbb_wind.dir,'+');
%         scatter(datenum(wind.time),wind.dir,'+');
scatter(swell.time,swell.dir,'go');
ylabel('Direction (°)');
ylim([0 360]); % weird outliers sometimes...

fg(5) = subplot(5,1,5);
yyaxis left
plot(swell.time,swell.hgt,'+');
ylabel('H_s (m)');
hold on
yyaxis right
scatter(swell.time,swell.per,'.');
ylabel('Period (s)');
ylim([0 20]);

linkaxes(fg,'x');
xlim([min(window_times_1) max(window_times_1)]);

figure
for ss = 1:n_sensors
    hold on
    eval(['wind_spd_mapped = interp1(datenum(wind.time),wind.spd,datenum(window_times_' num2str(ss) '));']);
    eval(['plot(window_times_' num2str(ss) ',(Hs_wind_' num2str(ss) '.^2)./(wind_spd_mapped.^2));']);
end
xlabel('Datetime');
ylabel('Wind Speed^2 / H_s Wind^2');
legend(labels);

figure
for ss = 1:n_sensors
    hold on
    eval(['plot(window_times_' num2str(ss) ',Hs_wind_' num2str(ss) ');']);
end
xlabel('Datetime');
ylabel('H_s from Wind');
legend(labels);

figure
for ss = 1:n_sensors
    hold on
    eval(['plot(window_times_' num2str(ss) ',Hs_swell_' num2str(ss) ');']);
end
xlabel('Datetime');
ylabel('H_s from Swell');
legend(labels);

figure
for ss = 1:n_sensors
    hold on
    eval(['plot(window_times_' num2str(ss) ',Hs_igw_' num2str(ss) ');']);
end
xlabel('Datetime');
ylabel('H_s from IGW');
legend(labels);


% dEF_LL_to_SB = energy_flux_2 - energy_flux_1;
% dEF_SB_to_PN = energy_flux_3 - energy_flux_2;
% dEF_PN_to_TP = energy_flux_4 - energy_flux_3;
% 
% figure
% semilogx(freq,dEF_LL_to_SB);
% hold on
% semilogx(freq,dEF_SB_to_PN);
% semilogx(freq,dEF_PN_to_TP);
% xlabel('frequency (Hz)');
% ylabel('Energy Flux Difference (m^3/s)');
% title(['Change in Energy Flux Going Into the Bay, ' datestr(start_time) ' to ' datestr(end_time)]);
% legend('LL to SB','SB to PN','PN to TP');
