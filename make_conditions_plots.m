% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last edited 1 October 2019

get_ocean_conditions
% load tomales_rbrs.mat
% Loads labels{}, rbr_depths{}, rbr_times{}
figure

f(1) = subplot(2,1,2);
yyaxis right
ylabel('Wind Speed (m/s)');
hold on
scatter(swell_time,windsp,'.');
ylim([0 20]);
hold on
yyaxis left
plot(swell_time,wvht,'+--');
ylabel('H_s (m)');

f(2) = subplot(2,1,1);
yyaxis left
plot(rbr_times{1},rbr_depths_adjusted{1});
% hold on
% plot(rbr_times{6},rbr_depths_adjusted{6});
ylabel('Depth at LL (m)');
yyaxis right
scatter(swell_time,dwpd,'.');
ylabel('Dominant Wave Period (s)');

% Making Tidal Index
% tmp = rbr_depths_adjusted{1} - mean(rbr_depths_adjusted{1});
% tmp = lowpass(tmp,0.0002,2);
% [max_peaks,max_loc] = findpeaks(tmp);
% max_peaks = max_peaks(2:end); max_loc = max_loc(2:end);
% [min_peaks,min_loc] = findpeaks(-tmp);
% min_peaks = -min_peaks;
% tmp = rbr_times{1};
% tide_switch_times = tmp(floor((max_loc+min_loc)./2));
% [spring_neap_index,sni_loc] = findpeaks(max_peaks-min_peaks,'MaxPeakWidth',10);
% spring_neap_index = spring_neap_index./max(spring_neap_index);
% spring_neap_index_times = tide_switch_times(sni_loc);
% plot(spring_neap_index_times,spring_neap_index,'r:','LineWidth',2);

% f(3) = subplot(3,1,3);
% yyaxis right
% scatter(swell_time,windsp,'.');
% ylabel('Wind Speed (m/s)');
% hold on
% yyaxis left
% plot(pressure_time,bpressure,'--');
% ylabel('Pressure (mbar)');

linkaxes(f,'x')
xlim([rbr_times{1}(1) rbr_times{1}(end)]);