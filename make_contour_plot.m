% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last Edited 17 July 2020

%% Load Control Parameters & Dates
run date_controls.m
run controls.m

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

%% Setup

% Cutoffs between types of waves (in seconds)
max_period_wind = 4;
max_period_swell = 25; % Maryam used 20s, we use 25
max_period_igw = 300; % Maryam used 300s, we use 250

% How many points do you want to use for the moving average of S? 
n_smooth = 5;

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
% max_varpreserv_power = 3.5*10^-3; % good for TB data
% max_varpreserv_power = 3*10^-3; % good for BB data

if end_time - start_time < 0.75
    use_t_tide = 0;
    fprintf('Time window less than 18 hours: will not use t_tide to de-tide as it gets imprecise.\n');
    fprintf('Also, will not use tidal depth binning as not enough of a study period.\n');
    spectra_choice = 0;
end
if instance_length > window_length, error('Ensemble Instance length must be smaller than the averaging window length.'), end
if window_length < ea_spacing, error('Right now, code requires windows be longer than ensemble spacing.'), end
extra = '';
if adjust_pressure
    extra = [extra 'with attenuation'];
end

T = seconds(rbr_times{1}(2)-rbr_times{1}(1)); % seconds
fs = 1/T; % Hz

depth_signal = rbr_depths_adjusted{sensor_choice};
% depth_signal = rbr_depths_raw{sensor_choice};

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
%     [~,d_out] = t_tide(trimmed_depth_signal,'start time',start,'latitude',38.17,'output','none');
    % Latitude given is for Tomales Point, but the above doesn't seem to
    % work. 
    eta = trimmed_depth_signal - d_out;
else
    eta = trimmed_depth_signal;
end

if pre_smooth
    eta = movmean(eta,6);
end

% Flatten any bottoming-out chunks ...?
if min(trimmed_depth_signal) <= 0, fprintf('!!! Depth signal bottoms out. Take note! Setting values to 0 !!!\n'), end    
trimmed_depth_signal(trimmed_depth_signal <= 0.01) = 0;

mab = cmab./100;
max_tide_range = max(trimmed_depth_signal) - min(trimmed_depth_signal); % in m

%% Spectral Machinery

[W,big_average_S] = f_ensemble_average_spectra(T,eta,instance_length*60*60,n_smooth);
freq = W./(2*pi); % Define Frequency! 

if adjust_pressure
    big_average_S = Sd_to_Ss(trimmed_depth_signal,mab(sensor_choice),W,big_average_S);
    % Boosts energy in the high-frequency bins...
end

measurements_per_window = window_length*60*60*fs;

% Decision: EA spacing dominates, we will clip timeseries at end if the
% combination of parameters doesn't "fit" perfectly into the timeseries. 

n_windows = floor((length(eta) - measurements_per_window)/(ea_spacing*60*60*fs));

matrixS = zeros(n_windows,length(big_average_S));
matrixSfreq = zeros(size(matrixS));
% Each row in matrix_S is a power spectra for a window. 
% The earlier time is at the top, the latter time at the bottom. 

window_depths = zeros(1,n_windows);
window_slope = zeros(1,n_windows);
window_Hs_wind = zeros(1,n_windows);
window_Hs_swell = zeros(1,n_windows);
window_Hs_igw = zeros(1,n_windows);
window_Hs_seiche = zeros(1,n_windows);
window_m0_wind = zeros(1,n_windows);
window_m0_swell = zeros(1,n_windows);
window_m0_igw = zeros(1,n_windows);
window_m0_seiche = zeros(1,n_windows);
window_times = [datetime(2020,1,1,1,1,1)]; % To initiate it
running_S = zeros(size(big_average_S));

switch spectra_choice
    case 1
        high_tide_running_S = zeros(1,length(big_average_S));
        n_high_tide_S = 0;
        med_tide_running_S = zeros(1,length(big_average_S));
        n_med_tide_S = 0;
        low_tide_running_S = zeros(1,length(big_average_S));
        n_low_tide_S = 0;
    case 2
        morning_running_S = zeros(1,length(big_average_S));
        n_morning_S = 0;
        afternoon_running_S = zeros(1,length(big_average_S));
        n_afternoon_S = 0;
        evening_running_S = zeros(1,length(big_average_S));
        n_evening_S = 0;
    case 3
        flooding_running_S = zeros(1,length(big_average_S));
        n_flooding_S = 0;
        ebbing_running_S = zeros(1,length(big_average_S));
        n_ebbing_S = 0;
        slack_running_S = zeros(1,length(big_average_S));
        n_slack_S = 0;
    case 4
        flooding_running_S = zeros(1,length(big_average_S));
        n_flooding_S = 0;
        ebbing_running_S = zeros(1,length(big_average_S));
        n_ebbing_S = 0;
        low_running_S = zeros(1,length(big_average_S));
        n_low_S = 0;
        high_running_S = zeros(1,length(big_average_S));
        n_high_S = 0;
end


%% Loop to find state in the tidal cycle; necessary to pre-calculate. 
for nn = 0:n_windows-1
    si = nn*ea_spacing*60*60*fs + 1;
    ei = nn*ea_spacing*60*60*fs + measurements_per_window;
    
    [fitcoeff,~,~] = polyfit(datenum(trimmed_times(si:ei)),trimmed_depth_signal(si:ei),1);
    window_slope(nn+1) = fitcoeff(1);
    
    window_depths(nn+1) = mean(trimmed_depth_signal(si:ei));
end
if spectra_choice ~= 0
    [~,peak_locs] = findpeaks(window_depths);
    [~,trough_locs] = findpeaks(-window_depths);
end


%% Loop for the main machinery. 
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
    g = 9.80665;
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
    
    switch spectra_choice
        case 1
            if mean(trimmed_depth_signal(si:ei)) > min(trimmed_depth_signal)+(2/3)*max_tide_range
                high_tide_running_S = high_tide_running_S + S;
                n_high_tide_S = n_high_tide_S + 1;
            elseif mean(trimmed_depth_signal(si:ei)) < min(trimmed_depth_signal)+(1/3)*max_tide_range
                low_tide_running_S = low_tide_running_S + S;
                n_low_tide_S = n_low_tide_S + 1;
            else
                med_tide_running_S = med_tide_running_S + S;
                n_med_tide_S = n_med_tide_S + 1;
            end
        case 2   
            if hour(window_times(nn+1)) <= 5
                morning_running_S = morning_running_S + S;
                n_morning_S = n_morning_S + 1;
            elseif hour(window_times(nn+1)) > 5 && hour(window_times(nn+1)) < 18
                afternoon_running_S = afternoon_running_S + S;
                n_afternoon_S = n_afternoon_S + 1;
            else
                evening_running_S = evening_running_S + S;
                n_evening_S = n_evening_S + 1;
            end
        case 3
            if window_slope(nn+1) > 0.1
                flooding_running_S = flooding_running_S + S;
                n_flooding_S = n_flooding_S + 1;
            elseif window_slope(nn+1) < -0.1
                ebbing_running_S = ebbing_running_S + S;
                n_ebbing_S = n_ebbing_S + 1;
            else
                slack_running_S = slack_running_S + S;
                n_slack_S = n_slack_S + 1;
            end
        case 4
            [~,nearest_peak] = min((peak_locs-nn).^2); nearest_peak = peak_locs(nearest_peak);
            [~,nearest_trough] = min((trough_locs-nn).^2); nearest_trough = trough_locs(nearest_trough);
            
            if window_slope(nn+1) > 0.07
                flooding_running_S = flooding_running_S + S;
                n_flooding_S = n_flooding_S + 1;
            elseif window_slope(nn+1) < -0.07
                ebbing_running_S = ebbing_running_S + S;
                n_ebbing_S = n_ebbing_S + 1;
            elseif window_slope(nn+1) <= 0.07 && window_slope(nn+1) >= -0.07 && abs(nearest_trough-nn) <= 3
                low_running_S = low_running_S + S;
                n_low_S = n_low_S + 1;
            elseif window_slope(nn+1) <= 0.07 && window_slope(nn+1) >= -0.07 && abs(nearest_peak-nn) <= 3
                high_running_S = high_running_S + S;
                n_high_S = n_high_S + 1;
            else
                fprintf('Not clear how to categorize window %d.\n',nn);
            end
    end
    
    matrixS(nn+1,:) = S; % movmean(S,5); % Already got smoothed in f_ensemble_average_spectra
    running_S = running_S + S;
    matrixSfreq(nn+1,:) = matrixS(nn+1,:).*freq;
    
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
    window_Hs_wind(nn+1) = 4*sqrt(m0_wind);
    window_Hs_swell(nn+1) = 4*sqrt(m0_swell);
    window_Hs_igw(nn+1) = 4*sqrt(m0_igw);
    window_Hs_seiche(nn+1) = 4*sqrt(m0_seiche);
end

running_S = running_S./n_windows;

switch spectra_choice
    case 1
        high_tide_running_S = high_tide_running_S./n_high_tide_S;
        med_tide_running_S = med_tide_running_S./n_med_tide_S;
        low_tide_running_S = low_tide_running_S./n_low_tide_S;
    case 2
        morning_running_S = morning_running_S./n_morning_S;
        afternoon_running_S = afternoon_running_S./n_afternoon_S;
        evening_running_S = evening_running_S./n_evening_S;
    case 3
        flooding_running_S = flooding_running_S./n_flooding_S;
        ebbing_running_S = ebbing_running_S./n_ebbing_S;
        slack_running_S = slack_running_S./n_slack_S;
    case 4
        flooding_running_S = flooding_running_S./n_flooding_S;
        ebbing_running_S = ebbing_running_S./n_ebbing_S;
        high_running_S = high_running_S./n_high_S;
        low_running_S = low_running_S./n_low_S;
end


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

Hs_wind = 4*sqrt(m0_wind);
Hs_swell = 4*sqrt(m0_swell);
Hs_igw = 4*sqrt(m0_igw);
Hs_seiche = 4*sqrt(m0_seiche);
Hs_total = 4*sqrt(m0_total);

% fprintf('Hs wind = %f, Hs swell = %f, Hs IGW = %f, Hs Seiche = %f\n',Hs_wind,Hs_swell,Hs_igw,Hs_seiche);
fprintf('Hs wind = %f, Hs swell = %f, Hs IGW = %f, Hs TOTAL = %f\n',Hs_wind,Hs_swell,Hs_igw,Hs_total);


% Pie Chart Shit
% figure
% if include_seiche
%     pie([m0_wind m0_swell m0_igw m0_seiche].*(10^7),{'Wind','Swell','IGW','Seiche'});
% else
%     pie([m0_wind m0_swell m0_igw].*(10^7),{'Wind','Swell','IGW'});
% end
% title(['Energy (m^2) Budget Partition for ' labels{sensor_choice}]);

% Bar Chart Shit
figure
if include_seiche
    percent_wind = m0_wind/(m0_wind+m0_swell+m0_igw+m0_seiche);
    percent_swell = m0_swell/(m0_wind+m0_swell+m0_igw+m0_seiche);
    percent_igw = m0_igw/(m0_wind+m0_swell+m0_igw+m0_seiche);
    percent_seiche = m0_seiche/(m0_wind+m0_swell+m0_igw+m0_seiche);
    bar_labels = categorical({'Wind','Swell','IGW','Seiche'});
    bar_labels = reordercats(bar_labels,{'Wind','Swell','IGW','Seiche'});
    b = bar(bar_labels,[m0_wind,m0_swell,m0_igw,m0_seiche]);
    xtips = b.XEndPoints;
    ytips = b.YEndPoints;
    text(xtips,ytips,{num2str(percent_wind,2),num2str(percent_swell,2),num2str(percent_igw,2),num2str(percent_seiche,2)},'HorizontalAlignment','center','VerticalAlignment','bottom');
else
    percent_wind = m0_wind/(m0_wind+m0_swell+m0_igw);
    percent_swell = m0_swell/(m0_wind+m0_swell+m0_igw);
    percent_igw = m0_igw/(m0_wind+m0_swell+m0_igw);
    bar_labels = categorical({'Wind','Swell','IGW'});
    bar_labels = reordercats(bar_labels,{'Wind','Swell','IGW'});
    b = bar(categorical({'Wind','Swell','IGW'}),[m0_wind,m0_swell,m0_igw]);
    xtips = b.XEndPoints;
    ytips = b.YEndPoints;
    text(xtips,ytips,{num2str(percent_wind,2),num2str(percent_swell,2),num2str(percent_igw,2)},'HorizontalAlignment','center','VerticalAlignment','bottom');
end
title(['Energy (m^2) Values for ' labels{sensor_choice}]);
ylim([0 4*10^-5]);

%% Conditions

% get_conditions
load conditions.mat
load tbb_wind.mat
load buoy_spectra.mat

% Spring-Neap Index
% Every m hours, compute index based off min-max...?
% Kinda janky to smooth it out enough but let's give it a shot. 
% m = 8;
% n_chunks = floor(length(trimmed_depth_signal)./(m*60*60*fs));
% diffs = zeros(1,n_chunks);
% sni.time = [datetime(2020,1,1,1,1,1)];
% for ii = 0:n_chunks-1
%     si = floor(ii*m*60*60*fs + 1);
%     ei = floor((ii+1)*m*60*60*fs);
%     diffs(ii+1) = max(trimmed_depth_signal(si:ei)) - min(trimmed_depth_signal(si:ei));
%     sni.time(ii+1) = trimmed_times(ei);
% end
% tmp_diffs = movmean(movmean(diffs,3),8);
% sni.value = (tmp_diffs-min(tmp_diffs))./(max(tmp_diffs)-min(tmp_diffs));

%% Energy Over Time & Correlations Stuff

if do_energy_over_time

    % Interpolate conditions onto energy times
    % Note, these don't account for direction of wind or swell. 
%     wind_for_corr = interp1(datenum(tbb_wind.time),tbb_wind.spd,datenum(ea_times)); % Note, using TBB wind here
%     wind_for_corr = interp1(datenum(wind.time),wind.spd,datenum(ea_times));
%     swell_for_corr = interp1(datenum(swell.time),swell.hgt,datenum(ea_times));
%     sni_for_corr = interp1(datenum(sni.time),sni.value,datenum(ea_times));
    
%     figure
%     hold on
%     plot(datenum(ea_times),wind_energy./max(wind_energy));
%     plot(datenum(ea_times),swell_energy./max(swell_energy));
%     plot(datenum(ea_times),igw_energy./max(igw_energy));
%     plot(datenum(ea_times),seiche_energy./max(seiche_energy));
%     plot(datenum(trimmed_times),eta./max(eta),'g:');
%     xlim([datenum(ea_times(1)) datenum(ea_times(end))]);
%     title(['Normalized Energy from Wind, Swell, IGW, Seiche Bands at ' labels{sensor_choice}]);
%     legend('Sea','Swell','IGW','Seiche','water surface');
    
    % For cross- and auto-correlation...
%     R = corrcoef(wind_for_corr,m0_wind); R = R(2)
%     corrcoef(sni_for_corr(4:1347),igw_energy(4:1347))

%     [acf,lags,~] = autocorr(igw_energy);
%     plot(lags.*ea_spacing,acf);

    figure
    sgtitle([labels{sensor_choice} ' H_s (via m_0), ' datestr(start_time) ' to ' datestr(end_time) '. ' extra]);
    if visualize_conditions
        ff(1) = subplot(4,2,[1,2,3,4]);
    end

    plot(window_times,window_Hs_wind);
    hold on
    plot(window_times,window_Hs_swell);
    plot(window_times,window_Hs_igw);
    if include_seiche
        plot(window_times,window_Hs_seiche);
    end
    yyaxis right
    plot(trimmed_times,trimmed_depth_signal,':');
    if include_seiche
        legend('Wind','Swell','IGW','Seiche','Depth Signal');
    else
        legend('Wind','Swell','IGW','Depth Signal');
    end
    
    eval(['Hs_wind_' num2str(sensor_choice) ' = window_Hs_wind;']);
    eval(['Hs_swell_' num2str(sensor_choice) ' = window_Hs_swell;']);
    eval(['Hs_igw_' num2str(sensor_choice) ' = window_Hs_igw;']);
    eval(['window_times_' num2str(sensor_choice) ' = window_times;']);
    
    if visualize_conditions
        ff(2) = subplot(4,2,[5,6]);
        yyaxis right

        scatter(tbb_wind.time,tbb_wind.spd,'.');
%         scatter(wind.time,wind.spd,'.');
        ylim([0 20]);
        ylabel('Wind Speed (m/s)');
        yyaxis left
        hold on
        scatter(tbb_wind.time,tbb_wind.dir,'+');
%         scatter(wind.time,wind.dir,'+');
        scatter(swell.time,swell.dir,'go');
        ylabel('Direction (°)');
        ylim([0 360]); % weird outliers sometimes...

        ff(3) = subplot(4,2,[7,8]);
        yyaxis left
        plot(swell.time,swell.hgt,'+');
        ylabel('H_s (m)');
        hold on
%         plot(datenum(sni.time),sni.value.*5,'g-');
        plot(trimmed_times,trimmed_depth_signal,'k-');
        yyaxis right
        scatter(swell.time,swell.per,'.');
        ylabel('Period (s)');
        ylim([0 20]);
        
        linkaxes(ff,'x');
        xlim([min(trimmed_times) max(trimmed_times)]);
    end
end

%% Plotting

logSt = logSt'; % Transpose! 
matrixSfreq = matrixSfreq';

if make_spectra_plot
    figure
    sgtitle([labels{sensor_choice} ' Spectra, ' num2str(ea_spacing) '/' num2str(window_length) '/' num2str(instance_length) '. ' extra]);

    if visualize_conditions
        ff(1) = subplot(4,2,[1,2,3,4]);
    end

    % It is good to start at 2 to avoid inf. 
    if variance_preserving
        % Impose high & low values to fix scale bar. 
%         matrixSfreq(end,end) = max_varpreserv_power;
%         matrixSfreq(end,end-1) = 0;
        contourf(datenum(window_times(:)),logfreq(n_leakage_ignore:end),matrixSfreq(n_leakage_ignore:end,:),15,'LineColor','none');
    else
        % Impose high & low values to fix scale bar. 
        logSt(2,end) = -1;
        logSt(2,end-1) = -9;
        contourf(datenum(window_times(:)),logfreq(n_leakage_ignore:end),logSt(n_leakage_ignore:end,:),15,'LineColor','none'); % 9 is pretty good, or 15
    end
    c = colorbar('east');
    if use_residual_spectra
        c.Label.String = 'Residual Power Density (m^2/Hz)';
    else
        c.Label.String = 'Power Density (m^2/Hz)';
    end

    % datetick('x',21);
    xlim([datenum(window_times(1)) datenum(window_times(end))]);
    % new_ticks = linspace(datenum(ea_times(1)),datenum(ea_times(end)),8);
    % xticks(new_ticks);
    % new_tick_labels = {};
    % for mm = 1:length(new_ticks)
    %     new_tick_labels{mm} = datestr(new_ticks(mm));
    % end
    % xticklabels(new_tick_labels);
    % xtickangle(25);
    ylabel('Log_{10} Frequency (Hz)'); 
    xlabel('Time (Center of Interval)');
    zlabel('Log Power Density (m^2/Hz)');

    if visualize_conditions
        ff(2) = subplot(4,2,[5,6]);
        yyaxis right

        scatter(datenum(tbb_wind.time),tbb_wind.spd,'.');
%         scatter(datenum(wind.time),wind.spd,'.');
        ylim([0 20]);
        ylabel('Wind Speed (m/s)');
        yyaxis left
        hold on
        scatter(datenum(tbb_wind.time),tbb_wind.dir,'+');
%         scatter(datenum(wind.time),wind.dir,'+');
        scatter(datenum(swell.time),swell.dir,'go');
        ylabel('Direction (°)');
        ylim([0 360]); % weird outliers sometimes...
    %     xlim([datetime(start_time,'ConvertFrom','datenum') datetime(end_time,'ConvertFrom','datenum')]);
        xlim([start_time end_time]);

        ff(3) = subplot(4,2,[7,8]);
        yyaxis left
        plot(datenum(swell.time),swell.hgt,'+');
        ylabel('H_s (m)');
        hold on
%         plot(datenum(sni.time),sni.value.*5,'g-');
        plot(datenum(trimmed_times),trimmed_depth_signal,'k-');
        yyaxis right
        scatter(datenum(swell.time),swell.per,'.');
        ylabel('Period (s)');
        ylim([0 20]);
    %     xlim([datetime(start_time,'ConvertFrom','datenum') datetime(end_time,'ConvertFrom','datenum')]);
        xlim([start_time end_time]);
        
        linkaxes(ff,'x');
        
        wave_period_mapped = interp1(datenum(swell.time),swell.per,datenum(window_times'));
        wave_height_mapped = interp1(datenum(swell.time),swell.hgt,datenum(window_times'));
        wave_direction_mapped = interp1(datenum(swell.time),double(swell.dir),datenum(window_times')).*(2*pi/360);
        
        wind_speed_mapped = interp1(datenum(tbb_wind.time),tbb_wind.spd,datenum(window_times'));
        wind_dir_mapped = interp1(datenum(tbb_wind.time),tbb_wind.dir,datenum(window_times'));
        wind_direction_NW = wave_direction_mapped > 4.89;
        eval(['Hs_wind_timeseries = Hs_wind_' num2str(sensor_choice) ';']);
        Hs_wind_timeseries = Hs_wind_timeseries';
        
        figure
        sgtitle([labels{sensor_choice} ' Conditions from ' num2str(ea_spacing) '/' num2str(window_length) '/' num2str(instance_length) '. ']);
        hh(1) = subplot(1,2,1);
        polarscatter(wave_direction_mapped,wave_height_mapped,30,window_m0_igw,'filled');
        hold on
        polarscatter(0,0,30,0);
        polarscatter(0,0,30,3*10^-5);
        colormap cool
        pax = gca;
        pax.ThetaDir = 'clockwise';
        pax.ThetaZeroLocation = 'top';
        c = colorbar;
        c.Label.String = 'm_0 from IGW';
        title('Dominant Wave Direction (theta) and Offshore Wave Height (radius)');
        
        hh(2) = subplot(1,2,2);
        polarscatter(wave_direction_mapped,wave_period_mapped,30,window_m0_igw,'filled');
        hold on
        polarscatter(0,0,30,0);
        polarscatter(0,0,30,3*10^-5);
        colormap cool
        pax = gca;
        pax.ThetaDir = 'clockwise';
        pax.ThetaZeroLocation = 'top';
        c = colorbar;
        c.Label.String = 'm_0 from IGW';
        title('Dominant Wave Direction (theta) and Offshore Wave Period (radius)');
        
        % Wind Speeds
        U = 0:0.1:20; % m/s
        U_A = 0.71*U.^1.23; % per SPM 3-28a
        F = 13; % km

%         H_m0_spm_fl = (U_A.^2).*((g^-1)*1.6*(10^-3)*(g*F./(U_A.^2)).^(1/2));
%         H_m0_spm_fd = (U_A.^2).*((g^-1)*2.433*10^-1);
%         T_m0_spm_fl = (U_A).*((g^-1)*2.857*(10^-1)*(g*F./(U_A.^2)).^(1/3));
        H_m0_spm_fl = (5.112*10^-4).*U_A.*(F^0.5);
        H_m0_spm_fd = (2.482*20^-2).*(U_A.^2);
        
        n_points_to_average = 4;
        
        figure
%         scatter(window_depths,wind_speed_mapped,30,Hs_wind_timeseries,'filled');
%         scatter(wind_speed_mapped,Hs_wind_timeseries,30,'filled');
        tmp1 = wind_speed_mapped.*wind_direction_NW;
        tmp2 = Hs_wind_timeseries.*wind_direction_NW;
        wind_speed_mapped_left_movmean = movmean(tmp1,[n_points_to_average,0]);
        scatter(tmp1, tmp2, 25, wind_speed_mapped_left_movmean,'filled');
        fprintf('In Wind plot, only plotting points with NWish winds.\n');
        colormap cool
        c = colorbar;
        c.Label.String = [num2str(n_points_to_average) '-Window Left-Handed Mean Wind Speed'];
        
        hold on
        plot(U,H_m0_spm_fl);
        plot(U,H_m0_spm_fd);

        xlabel('Wind Speed (m/s)');
        ylabel('H_s in Wind Chop (m)');
        ylim([0 0.3]);
        xlim([0 25]);
%         xlabel('Water Depth (m)');
        legend('Data','Fetch-Limited','Fully-Developed');
        title(['Wind Chop Development at ' labels{sensor_choice} ', ' num2str(instance_length) '-long instances.']);
    end

end

%% "Big Spectrum" Plotting

start_color = [149,111,49]./255; % light brown
end_color = [51,237,233]./255; % aqua
colors = [linspace(start_color(1),end_color(1),3)',linspace(start_color(2),end_color(2),3)',linspace(start_color(3),end_color(3),3)'];

if visualize_big_spectra
    figure
    if variance_preserving
        if use_median_power
            medianS = median(matrixS);
            semilogx(freq(n_leakage_ignore:end),medianS(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'k');
        else
            semilogx(freq(n_leakage_ignore:end),running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'k');
        end
        hold on
        
        switch spectra_choice
            case 1
                semilogx(freq(n_leakage_ignore:end),high_tide_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'Color',colors(1,:));
                semilogx(freq(n_leakage_ignore:end),med_tide_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'Color',colors(2,:));
                semilogx(freq(n_leakage_ignore:end),low_tide_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'Color',colors(3,:));
                legend('Total Running','High Tide','Med Tide','Low Tide');
            case 2
                semilogx(freq(n_leakage_ignore:end),morning_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'r');
                semilogx(freq(n_leakage_ignore:end),afternoon_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'g');
                semilogx(freq(n_leakage_ignore:end),evening_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'b');
                legend('Total Running','Morning','Afternoon','Evening');
            case 3
                semilogx(freq(n_leakage_ignore:end),flooding_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'r');
                semilogx(freq(n_leakage_ignore:end),ebbing_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'g');
                semilogx(freq(n_leakage_ignore:end),slack_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'b');
                legend('Total Running','Flooding','Ebbing','Slack');
            case 4
                semilogx(freq(n_leakage_ignore:end),flooding_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'r');
                semilogx(freq(n_leakage_ignore:end),ebbing_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'g');
                semilogx(freq(n_leakage_ignore:end),high_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'b');
                semilogx(freq(n_leakage_ignore:end),low_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'m');
                legend('Total Running','Flooding','Ebbing','High Slack','Low Slack');
            otherwise
                % Nothing! 
        end
        
        ylabel('Variance Per Hz');
%         tmp_max = max([high_tide_running_S med_tide_running_S low_tide_running_S]);
        axis([min(freq),max(freq),0,max_varpreserv_power]);
        % 2*10^-3 is a good alternative upper bound
        
    else
%         loglog(freq,big_average_S,'b');
%         hold on
        loglog(freq(n_leakage_ignore:end),running_S(n_leakage_ignore:end),'b');
        axis([min(freq),max(freq),10^-8,10^-2]);
        ylabel('Depth Power Density (m^2/Hz)');
    end
    xlabel('Frequency (Hz)');
    title(['Mean of ' num2str(instance_length) '-hr Spectra at ' labels{sensor_choice} ' from ' datestr(start_time) ' to ' datestr(end_time)]);
    
    % OK actually we're just gonna fit to the juicy middle part...
    [~,si] = min(abs(freq - 10^-2));
    [~,ei] = min(abs(freq - 10^-1));
    
    % Or maybe just 0.05 Hz to the end, per Hughes et al 2014
%     [~,si] = min(abs(freq - 0.05)); % to end

    fitcoeff = polyfit(logfreq(si:ei),logBigS(si:ei),1);
    
    fprintf('Spectra curve slope is %f\n',fitcoeff(1));

    
end

%% Energy Flux & Dissipation Stuff

run energy_flux.m

eval(['energy_flux_' num2str(sensor_choice) ' = running_S.*Cg;']);
eval(['running_S_' num2str(sensor_choice) ' = running_S;']);
% DON'T FORGET... 
% THIS IS ENERGY FLUX PER FREQUENCY, NOT YET FULLY DIMENSIONED EITHER

%% Tidal Depth Histogram

% figure
% histogram(window_depths,10,'Normalization','probability');
% xlabel('Tidal Depth of Window');
% ylabel('Proportion');
% title(labels{sensor_choice});

%% Closing

fprintf('\n\n');
