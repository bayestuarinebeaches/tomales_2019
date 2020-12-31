% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last Edited 14 Sept 2020

%% Load Control Parameters & Dates
run controls.m
run date_controls.m

%% Liftoff Announcements

fprintf(['*** RUNNING CODE FOR ' labels{sensor_choice} ' with %.2f, %.2f, %.2f. ***\n'],ea_spacing,window_length,instance_length);
if adjust_pressure
    fprintf('Accounting for pressure attenuation. ');
else
    fprintf('No pressure adjustment. ');
end
if use_t_tide
    fprintf('Using t_tide. ');
else
    fprintf('Dumb de-tiding. ');
end
if use_windowing
    fprintf('Windowing. ');
else
    fprintf('No windowing. ');
end
if low_pass_filter
    fprintf('Low-pass filtering. ');
end
fprintf('\n');
switch wind_option
    % 1 uses BOON, 2 uses TBB, 3 uses HIOC, 4 uses TBOC
    case 1
        fprintf('Using BOON Wind data.\n');
    case 2
        fprintf('Using TBB Wind data.\n');
    case 3
        fprintf('Using HIOC Wind data.\n');
    case 4
        fprintf('Using TBOC Wind data.\n');
    case 5
        fprintf('Using TBB wind speeds, BOON wind directions.\n');
end

%% Setup Parameters

% Cutoffs between types of waves (in seconds)
max_period_wind = 4;
max_period_swell = 25; % Maryam used 20s, we use 25
max_period_igw = 300; % Maryam used 300s

% How many points do you want to use for the moving average of S? 
n_smooth = 3;

% How many low-frequency bins will we ignore?
% Minimum needs to be 1 to avoid weird plotting of 0/inf
% Usually I'll use 3
n_leakage_ignore = 3;
% Then we step it up one for indexing...
n_leakage_ignore = n_leakage_ignore + 1;

% What frequency to dictate low-pass filtering? 
% low_pass_freq = 0.7; % Hz
% low_pass_freq = 0.08; % Trying to filter swell to evaluate tidal height in BB

max_varpreserv_power = 1*10^-4; % good for TB data
% max_varpreserv_power = 3.5*10^-3; % also good for TB data
% max_varpreserv_power = 3*10^-3; % good for BB data

g = 9.80665;

ws = (((D50s(sensor_choice)*10^-3)^2)/18)*(1650/(1.31*10^-3))*g; % Settling velocity using D50; note converting mm to m
% Note this comes from Dyer 1986 p109; assuming 1.31*10^-3 as mu

if end_time - start_time < 0.75
    use_t_tide = 0;
    fprintf('Time window less than 18 hours: will not use t_tide to de-tide as it gets imprecise.\n');
    fprintf('Also, tidal depth binning probably ineffective.\n');
end
if instance_length > window_length, error('Ensemble Instance length must be smaller than the averaging window length.'), end
if window_length < ea_spacing, error('Right now, code requires windows be longer than ensemble spacing.'), end
extra = '';
if adjust_pressure
    extra = [extra 'with attenuation'];
end

T = seconds(rbr_times{1}(2)-rbr_times{1}(1)); % seconds
fs = 1/T; % Hz

%% Set up the Signal
depth_signal = rbr_depths_adjusted{sensor_choice}; % depth_signal = rbr_depths_raw{sensor_choice} % for testing

times = datenum(rbr_times{sensor_choice});
start_index = find(times == start_time);
end_index = find(times == end_time);

if isempty(start_index) || isempty(end_index)
    error('Start or End Time index not found. Check dates and such.');
end
trimmed_depth_signal = depth_signal(start_index:end_index); 
trimmed_times = rbr_times{sensor_choice}(start_index:end_index); % Comes as datetimes

% Low-Pass Filtering
if low_pass_filter
    trimmed_depth_signal = lowpass(trimmed_depth_signal,low_pass_freq,fs);
end

% Remove the tidal signal using t_tide (https://www.eoas.ubc.ca/~rich/)
if use_t_tide
    [~,d_out] = t_tide(trimmed_depth_signal,T/3600,'output','none');
    % Latitude given is for Tomales Point, but the above doesn't seem to
    % work. 
    eta = trimmed_depth_signal - d_out;
else
    eta = trimmed_depth_signal;
end

if pre_smooth
    eta = movmean(eta,n_smooth);
end

% Flatten any bottoming-out chunks.
if min(trimmed_depth_signal) <= 0
    fprintf('!!! Depth signal bottoms out. Take note! Setting values to 0 !!!\n') 
    trimmed_depth_signal(trimmed_depth_signal <= 0.01) = 0;
end

mab = cmab./100;
max_tide_range = max(trimmed_depth_signal) - min(trimmed_depth_signal); % in m

%% Set Up Spectral Data

[W,big_average_S] = f_ensemble_average_spectra(T,eta,instance_length*60*60,n_smooth);
freq = W./(2*pi); % Defines Frequency! 

if adjust_pressure
    big_average_S = Sd_to_Ss(trimmed_depth_signal,mab(sensor_choice),W,big_average_S);
    % Boosts energy in the high-frequency bins at large water depths. 
end

measurements_per_window = window_length*60*60*fs;

% Decision: EA spacing dominates, we will clip timeseries at end if the
% combination of parameters doesn't "fit" perfectly into the timeseries. 

n_windows = floor((length(eta) - measurements_per_window)/(ea_spacing*60*60*fs));

matrixS = zeros(n_windows,length(big_average_S));
matrixSfreq = zeros(size(matrixS));
matrixSCg = zeros(size(matrixS));
% Each row in matrix_S is a power spectra for a window. 
% The earlier time is at the top, the latter time at the bottom. 

h = mean(trimmed_depth_signal) + cmab(sensor_choice)/100;
% Using average depth as sensor depth
% h = 11; % m at the mouth, approximately

kh = qkhfs(W,h); % Oy, taking in a lot of time here...
k = kh./h; % a bit sloppy
Cp = ((g./k).*tanh(kh)).^0.5;
Cg = Cp.*(0.5+kh./sinh(2.*kh));

c_f = 0.01; % Per Tucker & Pitt p321
rho = 1021; % arbitrary guess

window_depths = zeros(1,n_windows);
window_slope = zeros(1,n_windows);
window_Hs_wind = zeros(1,n_windows);
window_Hs_swell = zeros(1,n_windows);
window_Hs_igw = zeros(1,n_windows);
window_Hs_seiche = zeros(1,n_windows);
window_Hs_total = zeros(1,n_windows);
window_m0_wind = zeros(1,n_windows);
window_m0_swell = zeros(1,n_windows);
window_m0_igw = zeros(1,n_windows);
window_m0_seiche = zeros(1,n_windows);
window_m0_total = zeros(1,n_windows);
window_times = [datetime(2020,1,1,1,1,1)]; % To initiate it
window_Tp = zeros(1,n_windows);
window_m2_total = zeros(1,n_windows);
window_AvgT = zeros(1,n_windows); % "APD" on NOAA NDBC terminiology
window_nondimparameter = zeros(1,n_windows);
window_ub = zeros(1,n_windows);
window_taub = zeros(1,n_windows);
window_taustar = zeros(1,n_windows);

running_S = zeros(size(big_average_S));

flooding_running_S = zeros(1,length(big_average_S));
n_flooding_S = 0;
ebbing_running_S = zeros(1,length(big_average_S));
n_ebbing_S = 0;
low_running_S = zeros(1,length(big_average_S));
n_low_S = 0;
high_running_S = zeros(1,length(big_average_S));
n_high_S = 0;


%% Loop to find state in the tidal cycle; necessary to pre-calculate. 
for nn = 0:n_windows-1
    si = nn*ea_spacing*60*60*fs + 1;
    ei = nn*ea_spacing*60*60*fs + measurements_per_window;
    
    [fitcoeff,~,~] = polyfit(datenum(trimmed_times(si:ei)),trimmed_depth_signal(si:ei),1);
    window_slope(nn+1) = fitcoeff(1);
    
%     figure
%     plot(trimmed_times(si:ei),eta(si:ei));
%     title(['instance ' num2str(nn)]);
%     fprintf('figure %d has slope %f.\n',nn+1,window_slope(nn+1));
    
    window_depths(nn+1) = mean(trimmed_depth_signal(si:ei));
end

[~,peak_locs] = findpeaks(window_depths);
[~,trough_locs] = findpeaks(-window_depths);


%% Loop to find spectra of windows. 
for nn = 0:n_windows-1
    si = nn*ea_spacing*60*60*fs + 1;
    ei = nn*ea_spacing*60*60*fs + measurements_per_window;
    
    if ei > length(eta), error('Encountered an issue with indexing the signal.'), end
    
    ensemble = eta(si:ei);
    % Simple detrending if instances are short. 
    if instance_length < 1
        ensemble = detrend(ensemble,1); 
    elseif instance_length >= 1 && instance_length < 3
        ensemble = detrend(ensemble,2);
    end
        
    % Find index of high-frequency cutoff
    hfc = sqrt(g/(4*pi*(window_depths(nn+1))));
    
    if use_windowing
        % Uses Tukey Window, per recommendation here: 
        % https://arena-attachments.s3.amazonaws.com/5331399/0509f4edab72d8fb7e4628a7f6fcea3a.pdf?1571926476
        [W,S] = f_ensemble_average_spectra(T,ensemble.*tukeywin(length(ensemble)),instance_length*60*60,n_smooth);
    else
        [W,S] = f_ensemble_average_spectra(T,ensemble,instance_length*60*60,n_smooth);
    end
     
    window_times(nn+1) = trimmed_times(nn*ea_spacing*60*60*fs+1+floor(measurements_per_window/2)); % These times are at center of interval.
    
    if adjust_pressure
        S = Sd_to_Ss(ensemble,mab(sensor_choice),W,S);
    end
    
    [~,Tp_location] = max(freq.*S);
    window_Tp(nn+1) = 1/freq(Tp_location); % for Peak Period?note, this doesn't seem to work, needs more finessing. ltwp 21/12/2020
    
    [~,nearest_peak] = min((peak_locs-nn).^2); nearest_peak = peak_locs(nearest_peak);
    [~,nearest_trough] = min((trough_locs-nn).^2); nearest_trough = trough_locs(nearest_trough);

    if window_slope(nn+1) > window_slope_cutoff
        flooding_running_S = flooding_running_S + S;
        n_flooding_S = n_flooding_S + 1;
    elseif window_slope(nn+1) < -1*window_slope_cutoff
        ebbing_running_S = ebbing_running_S + S;
        n_ebbing_S = n_ebbing_S + 1;
    elseif abs(window_slope(nn+1)) <= window_slope_cutoff && abs(nearest_trough-nn) <= 3
        low_running_S = low_running_S + S;
        n_low_S = n_low_S + 1;
    elseif abs(window_slope(nn+1)) <= window_slope_cutoff && abs(nearest_peak-nn) <= 3
        high_running_S = high_running_S + S;
        n_high_S = n_high_S + 1;
    else
        fprintf('Not clear how to categorize window %d.\n',nn);
    end
    
    matrixS(nn+1,:) = S; % Already got smoothed in f_ensemble_average_spectra
    running_S = running_S + S;
    matrixSfreq(nn+1,:) = matrixS(nn+1,:).*freq;
    matrixSCg(nn+1,:) = matrixS(nn+1,:).*Cg;
    
    % Wind 
    [~,si] = min(abs(freq - 1/max_period_wind));
    [~,ei] = min(abs(freq - hfc));
    m0_wind = trapz(freq(si:ei),S(si:ei));

    % Swell
    [~,si] = min(abs(freq - 1/max_period_swell));
    [~,ei] = min(abs(freq - 1/max_period_wind));
    m0_swell = trapz(freq(si:ei),S(si:ei));

    % IGW 
    [~,si] = min(abs(freq - 1/max_period_igw)); 
    [~,ei] = min(abs(freq - 1/max_period_swell));
    m0_igw = trapz(freq(si:ei),S(si:ei));
    
    % Seiche
    [~,ei] = min(abs(freq - 1/max_period_igw));
    m0_seiche = trapz(freq(2:ei),S(2:ei));

    window_m0_wind(nn+1) = m0_wind;
    window_m0_swell(nn+1) = m0_swell;
    window_m0_igw(nn+1) = m0_igw;
    window_m0_seiche(nn+1) = m0_seiche;
    window_m0_total(nn+1) = m0_wind + m0_swell + m0_igw + m0_seiche;
    [~,ei] = min(abs(freq - hfc));
    window_m2_total(nn+1) = trapz(freq(2:ei).^2,S(2:ei));
    window_AvgT(nn+1) = sqrt(window_m0_total(nn+1)./window_m2_total(nn+1));
    window_Hs_total(nn+1) = 4*sqrt(window_m0_total(nn+1));
    window_Hs_wind(nn+1) = 4*sqrt(m0_wind);
    window_Hs_swell(nn+1) = 4*sqrt(m0_swell);
    window_Hs_igw(nn+1) = 4*sqrt(m0_igw);
    window_Hs_seiche(nn+1) = 4*sqrt(m0_seiche);
    window_nondimparameter(nn+1) = ws*window_AvgT(nn+1)/window_depths(nn+1);
    
    temp_kh = qkhfs(2*pi*1/window_AvgT(nn+1),window_depths(nn+1));
    window_ub(nn+1) = window_Hs_total(nn+1)*pi/(window_AvgT(nn+1)*sinh(temp_kh));
    window_taub(nn+1) = 0.5*0.025*rho*window_ub(nn+1)^2;
    window_taustar(nn+1) = window_taub(nn+1)/((2650-rho)*D50s(sensor_choice)*(10^-3)*g); % Shields Parameter
    
end

running_S = running_S./n_windows;
flooding_running_S = flooding_running_S./n_flooding_S;
ebbing_running_S = ebbing_running_S./n_ebbing_S;
high_running_S = high_running_S./n_high_S;
low_running_S = low_running_S./n_low_S;


%% Adjustment, Log Scale Mgmt
    
logfreq = log10(freq);
logS = log10(matrixS);
% logSfreq = log10(matrixSfreq);
logSt = zeros(size(logS));
logBigS = log10(big_average_S);
logRunningS = log10(running_S);
if use_residual_spectra
    for nn = 1:n_windows
        logSt(nn,:) = logS(nn,:) - logRunningS; 
        % Could subtract logBig_S instead, but it doesn't make as much sense with windowing and such. 
    end
else
    logSt = logS;
end

%% Conditions

% get_conditions
load conditions.mat
load buoy_spectra.mat
load hioc_wind.mat
load tboc_wind.mat

%% Overall Moments

% Wind 
[~,si] = min(abs(freq - 1/max_period_wind));
[~,ei] = min(abs(freq - hfc));
m0_wind = trapz(freq(si:ei),running_S(si:ei));

% Swell
[~,si] = min(abs(freq - 1/max_period_swell));
[~,ei] = min(abs(freq - 1/max_period_wind));
m0_swell = trapz(freq(si:ei),running_S(si:ei));

% IGW 
[~,si] = min(abs(freq - 1/max_period_igw)); 
[~,ei] = min(abs(freq - 1/max_period_swell));
m0_igw = trapz(freq(si:ei),running_S(si:ei));

% Seiche
[~,ei] = min(abs(freq - 1/max_period_igw));
m0_seiche = trapz(freq(2:ei),running_S(2:ei));

% total
[~,ei] = min(abs(freq - hfc));
m0_total = trapz(freq(n_leakage_ignore:ei),running_S(n_leakage_ignore:ei));
m2_total = trapz(freq(n_leakage_ignore:ei).^2,running_S(n_leakage_ignore:ei));

Hs_wind = 4*sqrt(m0_wind);
Hs_swell = 4*sqrt(m0_swell);
Hs_igw = 4*sqrt(m0_igw);
Hs_seiche = 4*sqrt(m0_seiche);
Hs_total = 4*sqrt(m0_total);
APD_total = sqrt(m0_total/m2_total);

% fprintf('Hs wind = %f, Hs swell = %f, Hs IGW = %f, Hs Seiche = %f\n',Hs_wind,Hs_swell,Hs_igw,Hs_seiche);
fprintf('Hs wind = %f, Hs swell = %f, Hs IGW = %f, Hs TOTAL = %f\n',Hs_wind,Hs_swell,Hs_igw,Hs_total);

%% Plot Generation

run plot_generation.m

%% Closing

fprintf('\n\n');
