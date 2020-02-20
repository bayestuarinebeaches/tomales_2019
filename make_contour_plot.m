% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last Edited 13 February 2020

%% Controls
%%%% ALL THREE OF THE NEXT VARIABLES ARE IN HOURS %%%%
ea_spacing = 3; % What is the spacing between each ensemble?
window_length = 24; % Each ensemble represents how many hours? 
instance_length = 1; % What is the length of each instance in an ensemble?

% 6/12/4 is a good combo
% 3/24/1 also useful

sensor_choice = 1;

cmab = [460, 20, 14, 47, 25, 116, 30];
    % Wall beach (5) and Tomales Beach (7) are guesses
    
% For FIRST Sensor Set
% start_time = datenum(2019,6,4,12,0,0);
% end_time = datenum(2019,7,16,0,0,0);
% load RBR_data/20190717/tomales_rbrs.mat

% For SECOND Sensor Set
start_time = datenum(2019,7,18,0,0,0);
end_time = datenum(2019,8,29,0,0,0);
load RBR_data/20190829/tomales_rbrs.mat

% For FOURTH Sensor Set
% start_time = datenum(2019,9,28,0,0,0);
% end_time = datenum(2019,11,24,0,0,0);
% load RBR_data/20191124/tomales_rbrs.mat


%% Setup

% tomales_rbrs has labels{}, rbr_depths{}, rbr_times{}, rbr_pressures{}, rbr_depths_adjusted{}

T = seconds(rbr_times{1}(2)-rbr_times{1}(1)); % s
fs = 1/T; % Hz

depth_signal = rbr_depths_adjusted{sensor_choice};
% depth_signal = rbr_depths{sensor_choice}; % Use if not yet adjusted
times = datenum(rbr_times{sensor_choice});
start_index = find(times == start_time);
end_index = find(times == end_time);
trimmed_depth_signal = depth_signal(start_index:end_index);
trimmed_times = rbr_times{sensor_choice}(start_index:end_index); % Comes as datetimes

% trimmed_depth_signal = pr_corr(trimmed_depth_signal,[],fs,cmab(sensor_choice));

%% Machinery

if instance_length > window_length, error('Ensemble Instance length must be smaller than the averaging window length.'), end
if window_length < ea_spacing, error('Right now, code requires windows be longer than ensemble spacing.'), end

measurements_per_window = window_length*60*60*fs;

% Decision: EA spacing dominates, we will clip timeseries at end if the
% combination of parameters doesn't "fit" perfectly. 

n_eas = ((end_time-start_time)*24/ea_spacing) - floor(window_length/ea_spacing);
[W,big_average_E] = f_ensemble_average_spectra(T,trimmed_depth_signal,instance_length*60*60,5);
matrix_E = zeros(n_eas,length(big_average_E));
% Each row in matrix_E is a power spectra for a place in time. 
% The earlier time is at the top, the latter time at the bottom. 

ea_times = [datetime(2020,1,1,1,1,1)]; % Just to initiate it

running_E = zeros(size(big_average_E));

for nn = 0:n_eas-1
    [W,E] = f_ensemble_average_spectra(T,trimmed_depth_signal((nn*ea_spacing*60*60*fs+1):(nn*ea_spacing*60*60*fs+measurements_per_window)),instance_length*60*60,5);
    ea_times(nn+1) = trimmed_times(nn*ea_spacing*60*60*fs+1+floor(measurements_per_window/2));
    freq = W./(2*pi);
    matrix_E(nn+1,:) = movmean(E,5);
    running_E = running_E + E;
end

running_E = running_E./n_eas;

logfreq = log10(freq);
logE = log10(matrix_E);
logEt = log10(matrix_E); 
logBigE = log10(big_average_E);
logRunning_E = log10(running_E);

%% Noise Removal
    
for nn = 1:n_eas-2
    logEt(nn,:) = logE(nn,:) - logBigE; % Could subtract logRunning_E instead
end

%% Conditions

get_ocean_conditions

%% Plotting

logEt = logEt'; % Transpose! 

figure
sgtitle([num2str(ea_spacing) '-hour spacing of ' labels{sensor_choice} ' spectra with ' num2str(window_length) '-hour Intervals with ' num2str(instance_length) '-hour Ensemble Windows']);

subplot(4,2,[1,2,3,4])
contourf(datenum(ea_times),logfreq(2:end),logEt(2:end,:),15,'LineColor','none'); % 9 is pretty good
c = colorbar('east');
c.Label.String = 'Residual Power Density (m^2/Hz)';

% OLD TICK MODIFICATION
% tmp = 0:n_eas/15:n_eas-1;
% xticks(tmp);
% tmp = tmp.*ea_spacing + window_length/2;
% tmp = start_time + tmp./24;
% xticklabels(datestr(tmp));

datetick('x',21);
xlim([datenum(ea_times(1)) datenum(ea_times(end))]);
ylabel('Log_{10} Frequency (Hz)'); 
xlabel('Time (Center of Window)');
zlabel('Log Power Density (m^2/Hz)');

subplot(4,2,[5,6]);
yyaxis right

scatter(swell_time,windsp,'.');
ylim([0 20]);
ylabel('Wind Speed (m/s)');
yyaxis left
hold on
scatter(swell_time,winddir,'o');
scatter(swell_time,wdir,'g+');
ylabel('Wind Direction (°)');
ylim([0 360]); % weird outliers sometimes...
xlim([datetime(start_time,'ConvertFrom','datenum') datetime(end_time,'ConvertFrom','datenum')]);

subplot(4,2,[7,8]);
yyaxis left
plot(swell_time,wvht,'+');
ylabel('H_s (m)');
hold on
plot(rbr_times{sensor_choice},depth_signal,'g-');
% plot(rbr_times{6},rbr_depths_adjusted{6});
% ylabel('Depth at LL (m)');
yyaxis right
scatter(swell_time,dwpd,'.');
ylabel('Dominant Wave Period (s)');
xlim([datetime(start_time,'ConvertFrom','datenum') datetime(end_time,'ConvertFrom','datenum')]);

% Quiver plot for winds & swell ? 




% figure
% loglog(freq,big_average_E,'b');
% hold on
% % loglog(freq,logrunning_E,'r');
% xlabel('Frequency (Hz)');
% ylabel('Power Density (m^2/Hz)');
% 
% title(['Arithmetic Mean of ' num2str(instance_length) '-hour Instance Spectra at ' labels{sensor_choice} ' from ' datestr(start_time) ' to ' datestr(end_time)]);


