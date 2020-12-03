% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last Edited 14 Sept 2020

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

%% Setup Parameters & Signal

% Cutoffs between types of waves (in seconds)
max_period_wind = 4;
max_period_swell = 25; % Maryam used 20s, we use 25
max_period_igw = 300; % Maryam used 300s, we use 250

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
window_m0_wind = zeros(1,n_windows);
window_m0_swell = zeros(1,n_windows);
window_m0_igw = zeros(1,n_windows);
window_m0_seiche = zeros(1,n_windows);
window_times = [datetime(2020,1,1,1,1,1)]; % To initiate it
window_Tp = zeros(1,n_windows);
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
    window_Tp(nn+1) = 1/freq(Tp_location);
    
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
    window_Hs_wind(nn+1) = 4*sqrt(m0_wind);
    window_Hs_swell(nn+1) = 4*sqrt(m0_swell);
    window_Hs_igw(nn+1) = 4*sqrt(m0_igw);
    window_Hs_seiche(nn+1) = 4*sqrt(m0_seiche);
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

if make_energy_bar_chart
    figure
    if include_seiche
        percent_wind = m0_wind/(m0_wind+m0_swell+m0_igw+m0_seiche);
        percent_swell = m0_swell/(m0_wind+m0_swell+m0_igw+m0_seiche);
        percent_igw = m0_igw/(m0_wind+m0_swell+m0_igw+m0_seiche);
        percent_seiche = m0_seiche/(m0_wind+m0_swell+m0_igw+m0_seiche);
        bar_labels = categorical({'Wind','Swell','IGW','Seiche'});
        bar_labels = reordercats(bar_labels,{'Wind','Swell','IGW','Seiche'});
        b = bar(bar_labels,[Hs_wind,Hs_swell,Hs_igw,Hs_seiche]);
        xtips = b.XEndPoints;
        ytips = b.YEndPoints;
        text(xtips,ytips,{num2str(percent_wind,2),num2str(percent_swell,2),num2str(percent_igw,2),num2str(percent_seiche,2)},'HorizontalAlignment','center','VerticalAlignment','bottom');
    else
        percent_wind = m0_wind/(m0_wind+m0_swell+m0_igw);
        percent_swell = m0_swell/(m0_wind+m0_swell+m0_igw);
        percent_igw = m0_igw/(m0_wind+m0_swell+m0_igw);
        bar_labels = categorical({'Wind','Swell','IGW'});
        bar_labels = reordercats(bar_labels,{'Wind','Swell','IGW'});
        b = bar(categorical({'Wind','Swell','IGW'}),[Hs_wind,Hs_swell,Hs_igw]);
        xtips = b.XEndPoints;
        ytips = b.YEndPoints;
%         text(xtips,ytips,{num2str(percent_wind,2),num2str(percent_swell,2),num2str(percent_igw,2)},'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',35);
        text(xtips,ytips,{num2str(Hs_wind,2),num2str(Hs_swell,2),num2str(Hs_igw,2)},'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
    title(['Energy Proportions at ' labels{sensor_choice}]);
    ylim([0 0.15]);
    ylabel('Wave Heights (m)');
end

%% Conditions

% get_conditions
load conditions.mat
load buoy_spectra.mat
load hioc_wind.mat
load tboc_wind.mat

%% Energy Over Time & Correlations Stuff

    % Interpolate conditions onto energy times
    % Note, these don't account for direction of wind or swell. 
%     wind_for_corr = interp1(datenum(tbb_wind.time),tbb_wind.spd,datenum(ea_times)); % Note, using TBB wind here
%     wind_for_corr = interp1(datenum(bml_wind.time),bml_wind.spd,datenum(ea_times));
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

if make_wave_height_plot
    figure
    sgtitle([labels{sensor_choice} ' H_s (via m_0), ' datestr(start_time) ' to ' datestr(end_time) '. ' extra]);

    % Curves of wave height over time. 
    ff(1) = subplot(4,2,[1,2,3,4]);

    plot(window_times,window_Hs_wind);
    hold on
    plot(window_times,window_Hs_swell);
    plot(window_times,window_Hs_igw);
    if include_seiche
        plot(window_times,window_Hs_seiche);
    end
%     yyaxis right
%     plot(trimmed_times,trimmed_depth_signal,':');
    if include_seiche
        legend('Wind','Swell','IGW','Seiche','Depth Signal');
    else
        legend('Wind','Swell','IGW','Depth Signal');
    end

    eval(['Hs_wind_' num2str(sensor_choice) ' = window_Hs_wind;']);
    eval(['Hs_swell_' num2str(sensor_choice) ' = window_Hs_swell;']);
    eval(['Hs_igw_' num2str(sensor_choice) ' = window_Hs_igw;']);
    eval(['window_times_' num2str(sensor_choice) ' = window_times;']);

    ff(2) = subplot(4,2,[5,6]);
    yyaxis right

    switch wind_option
        case 1
            scatter(bml_wind.spd_time,bml_wind.spd,'.');
        case 2
            scatter(tbb_wind.time,tbb_wind.spd,'.');
        case 3
            scatter(hioc_wind.time,hioc_wind.spd,'.');
        case 4
            scatter(tboc_wind.time,tboc_wind.spd,'.');
        case 5
            scatter(tbb_wind.time,tbb_wind.spd,'.');
    end
    ylim([0 20]);
    ylabel('Wind Speed (m/s)');
    yyaxis left
    hold on
    switch wind_option
        case 1
            scatter(bml_wind.dir_time,bml_wind.dir,'+');
        case 2
            scatter(tbb_wind.time,tbb_wind.dir,'+');
        case 3
            scatter(hioc_wind.time,hioc_wind.dir,'+');
        case 4
            scatter(tboc_wind.time,tboc_wind.dir,'+');
        case 5
            scatter(bml_wind.dir_time,bml_wind.dir,'+');
    end       
%     scatter(swell.time,swell.dir,'go');
    ylabel('Direction (°)');
    ylim([0 360]); % weird outliers sometimes...

    ff(3) = subplot(4,2,[7,8]);
    yyaxis left
    plot(swell.time,swell.hgt,'+');
    ylabel('H_s (m)');
    hold on
    plot(trimmed_times,trimmed_depth_signal,'k-');
    yyaxis right
    scatter(swell.time,swell.per,'.');
    ylabel('Period (s)');
    ylim([0 20]);

    linkaxes(ff,'x');
%     xlim([min(trimmed_times) max(trimmed_times)]);
end

%% Contour Plotting

logSt = logSt'; % Transpose! 
matrixSfreq = matrixSfreq';

if make_contour_graph
    figure
    sgtitle([labels{sensor_choice} ' Spectra, ' num2str(ea_spacing) '/' num2str(window_length) '/' num2str(instance_length) '. ' extra]);

    ff(1) = subplot(4,2,[1,2,3,4]);

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
        c.Label.String = 'Residual Energy Density (m^2/Hz)';
    else
        c.Label.String = 'Energy Density (m^2/Hz)';
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
    zlabel('Log Energy Density (m^2/Hz)');

    ff(2) = subplot(4,2,[5,6]);
    yyaxis right
    switch wind_option
        case 1
            scatter(datenum(bml_wind.spd_time),bml_wind.spd,'.');
        case 2
            scatter(datenum(tbb_wind.time),tbb_wind.spd,'.');
        case 3
            scatter(datenum(hioc_wind.time),hioc_wind.spd,'.');
        case 4
            scatter(datenum(tboc_wind.time),tboc_wind.spd,'.');
        case 5
            scatter(datenum(tbb_wind.time),tbb_wind.spd,'.');
    end
    
    ylim([0 20]);
    ylabel('Wind Speed (m/s)');
    yyaxis left
    hold on
    switch wind_option
        case 1
            scatter(datenum(bml_wind.dir_time),bml_wind.dir,'+');
        case 2
            scatter(datenum(tbb_wind.time),tbb_wind.dir,'+');
        case 3
            scatter(datenum(hioc_wind.time),hioc_wind.dir,'+');
        case 4
            scatter(datenum(tboc_wind.time),tboc_wind.dir,'+');
        case 5
            scatter(datenum(bml_wind.dir_time),bml_wind.dir,'+');
    end   
    scatter(datenum(swell.time),swell.dir,'go');
    ylabel('Direction (°)');
    ylim([0 360]); % weird outliers sometimes...
    xlim([start_time end_time]);

    ff(3) = subplot(4,2,[7,8]);
    yyaxis left
    plot(datenum(swell.time),swell.hgt,'+');
    ylabel('H_s (m)');
    hold on
    plot(datenum(trimmed_times),trimmed_depth_signal,'k-');
    yyaxis right
    scatter(datenum(swell.time),swell.per,'.');
    ylabel('Period (s)');
    ylim([0 20]);
    xlim([start_time end_time]);

    linkaxes(ff,'x');
    
end

wave_direction_mapped = interp1(datenum(swell.time),double(swell.dir),datenum(window_times'));
wave_period_mapped = interp1(datenum(swell.time),swell.per,datenum(window_times'));
wave_height_mapped = interp1(datenum(swell.time),swell.hgt,datenum(window_times'));

switch wind_option
    case 1
        wind_speed_mapped = interp1(datenum(bml_wind.spd_time),bml_wind.spd,datenum(window_times'));
        wind_dir_mapped = interp1(datenum(bml_wind.dir_time),double(bml_wind.dir),datenum(window_times'));
    case 2
        wind_speed_mapped = interp1(datenum(tbb_wind.time),tbb_wind.spd,datenum(window_times'));
        wind_dir_mapped = interp1(datenum(tbb_wind.time),tbb_wind.dir,datenum(window_times'));
    case 3
        wind_speed_mapped = interp1(datenum(hioc_wind.time),hioc_wind.spd,datenum(window_times'));
        wind_dir_mapped = interp1(datenum(hioc_wind.time),hioc_wind.dir,datenum(window_times'));
    case 4
        wind_speed_mapped = interp1(datenum(tboc_wind.time),tboc_wind.spd,datenum(window_times'));
        wind_dir_mapped = interp1(datenum(tboc_wind.time),tboc_wind.dir,datenum(window_times'));
    case 5
        % WEIRD UNIQUENESS ISSUE WITH BML WINDS...?
        [~,unique_index] = unique(bml_wind.dir_time);
        wind_speed_mapped = interp1(datenum(tbb_wind.time),tbb_wind.spd,datenum(window_times'));
        wind_dir_mapped = interp1(datenum(bml_wind.dir_time(unique_index)),bml_wind.dir(unique_index),datenum(window_times'));
end
wind_direction_NW = wind_dir_mapped > 280; % in degrees, 4.89 for radians
eval(['Hs_wind_timeseries = Hs_wind_' num2str(sensor_choice) ';']);
eval(['Hs_swell_timeseries = Hs_swell_' num2str(sensor_choice) ';']);
eval(['Hs_igw_timeseries = Hs_igw_' num2str(sensor_choice) ';']);
Hs_wind_timeseries = Hs_wind_timeseries';
Hs_swell_timeseries = Hs_swell_timeseries';
Hs_igw_timeseries = Hs_igw_timeseries';

if make_polar_plot
    figure
    sgtitle([labels{sensor_choice} ' Conditions from ' num2str(ea_spacing) '/' num2str(window_length) '/' num2str(instance_length) '. ']);
    hh(1) = subplot(1,2,1);
    polarscatter(wave_direction_mapped.*(2*pi/360),wave_height_mapped,25,window_Hs_igw,'filled');
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
    polarscatter(wave_direction_mapped.*(2*pi/360),wave_period_mapped,30,window_Hs_igw,'filled');
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
end

%% Wind Development Plotting

% Wind Speeds
U = 0:0.1:20; % m/s
U_A = 0.71*U.^1.23; % per SPM 3-28a
F = fetch(sensor_choice)*1000; % now in m

H_m0_spm_fl = (1.6*10^-3)*(g^-1).*(U_A.^2).*(g*F./(U_A.^2)).^(0.5);
T_m0_spm_fl = (2.857*10^-1)*(g^-1).*(U_A).*(g*F./(U_A.^2)).^(1/3);
t_m0_spm_fl = (6.88*10^1)*(g^-1).*(U_A).*(g*F./(U_A.^2)).^(2/3);
H_m0_spm_fd = (2.433*10^-1)*(g^-1).*(U_A.^2);
T_m0_spm_fd = 8.134*(g^-1).*(U_A);
t_m0_spm_fd = (7.15*10^4)*(g^-1).*(U_A);
% THE ABOVE NEEDS TO BE DOUBLE-CHECKED!!!

if make_chop_plot
    n_points_to_average = 4;
    
    tmp1 = wind_speed_mapped.*wind_direction_NW; % Good to try squared, too
    tmp2 = Hs_wind_timeseries.*wind_direction_NW;
    wind_speed_mapped_left_movmean = movmean(wind_speed_mapped,[n_points_to_average,0]);
    
    figure
    scatter(tmp1, tmp2, 25, window_depths,'filled');
    fprintf('In Chop plot, only plotting points with NWish winds.\n');
    hold on
    plot(U,H_m0_spm_fl);
    plot(U,H_m0_spm_fd);
    
    H_m0_spm_fl_attenuated = H_m0_spm_fl.*exp(- (2* 10^-3) * 1000); % k & x (m) approximation
    plot(U,H_m0_spm_fl_attenuated);
    
    colormap cool
    c = colorbar;
%     c.Label.String = [num2str(n_points_to_average) '-Window Left-Handed Mean Wind Speed'];
    c.Label.String = 'Average Water Depth for Window';
    xlabel('Wind Speed (m/s)');
    ylabel('H_s in Wind Chop (m)');
%     ylim([0 0.3]);
    xlim([0 25]);
%     set(gca,'yscale','log')
    legend('Data','Fetch-Limited','Fully-Developed','Attenuated Fetch-Limited');
    title(['Wind Chop Development at ' labels{sensor_choice} ', ' num2str(instance_length) '-long instances.']);
    
    figure
    scatter(tmp1,window_depths,35,tmp2,'filled');
    colormap cool
    c = colorbar;
    xlabel('Wind Speed (m/s)');
    ylabel('Water Depth (m)');
    c.Label.String = 'H_s Wind (m)';
    title(['Wind Chop Height at ' labels{sensor_choice} ', ' num2str(instance_length) '-long instances.']);
    
    figure
    scatter(window_depths,tmp2);
    xlabel('Window Water Depth (m)');
    ylabel('H_s in Wind Chop (m)');
    title(['Wind Chop Height at ' labels{sensor_choice} ', ' num2str(instance_length) '-long instances.']);
end

if make_swell_plot
    figure
    scatter(wave_height_mapped,window_depths,35,Hs_swell_timeseries,'filled');
    colormap cool
    c = colorbar;
    xlabel('Offshore Wave Height (m)');
    ylabel('Water Depth (m)');
    c.Label.String = 'H_s Swell (m)';
    title(['Swell Wave Height at ' labels{sensor_choice} ', ' num2str(instance_length) '-long instances.']);
    
    figure
    scatter(wave_height_mapped,window_depths,35,Hs_igw_timeseries,'filled');
    colormap cool
    c = colorbar;
    xlabel('Offshore Wave Height (m)');
    ylabel('Water Depth (m)');
    c.Label.String = 'H_s IGW (m)';
    title(['IGW Wave Height at ' labels{sensor_choice} ', ' num2str(instance_length) '-long instances.']);
end

if make_wind_start_plot
    figure
    switch wind_option
        case 1
            scatter(bml_wind.spd_time,bml_wind.spd,'.');
        case 2
            scatter(tbb_wind.time,tbb_wind.spd,'.');
        case 3
            scatter(hioc_wind.time, hioc_wind.spd,'.');
        case 4
            scatter(tboc_wind.time, tboc_wind.spd,'.');
        case 5
            scatter(tbb_wind.time,tbb_wind.spd,'.');
    end
    hold on
    ylim([0 25]);
    ylabel('Wind Speed (m/s)');
    yyaxis right
    plot(window_times,window_Hs_wind,'r-');
    ylabel('H_s due to Wind');
    xlim([min(window_times),max(window_times)]);
    xlabel('Date');
    title(['Wind Chop Development at ' labels{sensor_choice} ', ' num2str(instance_length) '-long instances.']);
    
end

if make_binned_chart
    % First Winds
    
    buckets = zeros(10,10);
    bucket_counts = zeros(10,10);
    [depths_binned,depths_edges] = discretize(window_depths,10);
    [windspeed_binned,windspeed_edges] = discretize(wind_speed_mapped,10);
    for mm = 1:n_windows
        buckets(windspeed_binned(mm),depths_binned(mm)) = buckets(windspeed_binned(mm),depths_binned(mm)) + window_m0_wind(mm);
        bucket_counts(windspeed_binned(mm),depths_binned(mm)) = bucket_counts(windspeed_binned(mm),depths_binned(mm)) + 1;
    end
    buckets = buckets./bucket_counts;
    buckets(isnan(buckets)) == 0;
    
    figure
    contourf(buckets,'LineColor','None')
    c = colorbar('east');
    xlabel('Binned Wind Speed');
    ylabel('Binned Water Depth');
    title(['Binning Wind and Depth at ' labels{sensor_choice} ', ' num2str(instance_length) '-long instances.']);
    
    % Now Swell
    
    buckets = zeros(10,10);
    bucket_counts = zeros(10,10);
    [swell_binned,windspeed_edges] = discretize(wave_height_mapped,10);
    for mm = 1:n_windows
        buckets(swell_binned(mm),depths_binned(mm)) = buckets(swell_binned(mm),depths_binned(mm)) + window_m0_swell(mm);
        bucket_counts(swell_binned(mm),depths_binned(mm)) = bucket_counts(swell_binned(mm),depths_binned(mm)) + 1;
    end
    buckets = buckets./n_windows;
    buckets(isnan(buckets)) == 0;
    
    figure
    contourf(buckets,'LineColor','None')
    c = colorbar('east');
    xlabel('Binned Offshore Wave Height');
    ylabel('Binned Water Depth');
    title(['Binning Swell and Depth at ' labels{sensor_choice} ', ' num2str(instance_length) '-long instances.']);
    
    
end

%% "Big Spectrum" Plotting

start_color = [149,111,49]./255; % light brown
end_color = [51,237,233]./255; % aqua
colors = [linspace(start_color(1),end_color(1),3)',linspace(start_color(2),end_color(2),3)',linspace(start_color(3),end_color(3),3)'];

if make_spectral_plot
    figure
    if variance_preserving
        if use_median_power
            medianS = median(matrixS);
            semilogx(freq(n_leakage_ignore:end),medianS(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'k');
        else
            semilogx(freq(n_leakage_ignore:end),running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'k');
        end
        hold on

        if ~plot_only_mean_spectrum
            semilogx(freq(n_leakage_ignore:end),flooding_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'r');
            semilogx(freq(n_leakage_ignore:end),ebbing_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'g');
            semilogx(freq(n_leakage_ignore:end),high_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'b');
            semilogx(freq(n_leakage_ignore:end),low_running_S(n_leakage_ignore:end).*freq(n_leakage_ignore:end),'m');
            legend('Total Running','Flooding','Ebbing','High Slack','Low Slack');
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
        ylabel('Depth Energy Density (m^2/Hz)');
    end
    xlabel('Frequency (Hz)');
    title(['Mean of ' num2str(instance_length) '-hr Spectra at ' labels{sensor_choice} ' from ' datestr(start_time) ' to ' datestr(end_time)]);

    % Bounds to fit... 
    [~,si] = min(abs(freq - 0.4)); % -2
    [~,ei] = min(abs(freq - hfc)); % -1 to  fit juicy middle park

    % Or maybe just 0.05 Hz to the end, per Hughes et al 2014
    %     [~,si] = min(abs(freq - 0.05)); % to end

    fitcoeff = polyfit(logfreq(si:ei),logRunningS(si:ei),1);

    fprintf('High-frequency spectral curve slope is %f\n',fitcoeff(1));
end


%% Energy Flux & Dissipation Stuff

run energy_flux.m

eval(['energy_flux_' num2str(sensor_choice) ' = running_S.*Cg;']);
eval(['running_S_' num2str(sensor_choice) ' = running_S;']);
eval(['flooding_running_S_' num2str(sensor_choice) ' = flooding_running_S;']);
eval(['ebbing_running_S_' num2str(sensor_choice) ' = ebbing_running_S;']);
eval(['high_running_S_' num2str(sensor_choice) ' = high_running_S;']);
eval(['low_running_S_' num2str(sensor_choice) ' = low_running_S;']);

if make_energy_flux_contour_graph
    figure
    sgtitle([labels{sensor_choice} ' Energy Flux, ' num2str(ea_spacing) '/' num2str(window_length) '/' num2str(instance_length) '. ' extra]);

    ff(1) = subplot(4,2,[1,2,3,4]);

    contourf(datenum(window_times(:)),logfreq(n_leakage_ignore:end),log(matrixSCg(:,n_leakage_ignore:end)'),15,'LineColor','none');

    c = colorbar('east');
    if use_residual_spectra
        c.Label.String = 'Residual Energy Density (m^2/Hz)';
    else
        c.Label.String = 'Ln Energy Flux (not yet dimensionalized)'; %(m^2/Hz)';
    end

    xlim([datenum(window_times(1)) datenum(window_times(end))]);    
    ylabel('Log_{10} Frequency (Hz)'); 
    xlabel('Time (Center of Interval)');

    ff(2) = subplot(4,2,[5,6]);
    yyaxis right
    switch wind_option
        case 1
            scatter(datenum(bml_wind.spd_time),bml_wind.spd,'.');
        case 2
            scatter(datenum(tbb_wind.time),tbb_wind.spd,'.');
        case 3
            scatter(datenum(hioc_wind.time),hioc_wind.spd,'.');
        case 4
            scatter(datenum(tboc_wind.time),tboc_wind.spd,'.');
        case 5
            scatter(datenum(tbb_wind.time),tbb_wind.spd,'.');
    end
    
    ylim([0 20]);
    ylabel('Wind Speed (m/s)');
    yyaxis left
    hold on
    switch wind_option
        case 1
            scatter(datenum(bml_wind.dir_time),bml_wind.dir,'+');
        case 2
            scatter(datenum(tbb_wind.time),tbb_wind.dir,'+');
        case 3
            scatter(datenum(hioc_wind.time),hioc_wind.dir,'+');
        case 4
            scatter(datenum(tboc_wind.time),tboc_wind.dir,'+');
        case 5
            scatter(datenum(bml_wind.dir_time),bml_wind.dir,'+');
    end   
    scatter(datenum(swell.time),swell.dir,'go');
    ylabel('Direction (°)');
    ylim([0 360]); % weird outliers sometimes...
    xlim([start_time end_time]);

    ff(3) = subplot(4,2,[7,8]);
    yyaxis left
    plot(datenum(swell.time),swell.hgt,'+');
    ylabel('H_s (m)');
    hold on
    plot(datenum(trimmed_times),trimmed_depth_signal,'k-');
    yyaxis right
    scatter(datenum(swell.time),swell.per,'.');
    ylabel('Period (s)');
    ylim([0 20]);
    xlim([start_time end_time]);

    linkaxes(ff,'x');
    
end

% DON'T FORGET... 
% THIS IS ENERGY FLUX PER FREQUENCY, NOT YET FULLY DIMENSIONED EITHER 
% (need to multiply by 0.5 * rho * g)

%% Closing

fprintf('\n\n');
