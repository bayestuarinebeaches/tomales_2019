% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last Edited 9 June 2020

%% Sensor Data

% Bring in labels{}, rbr_depths{}, rbr_times{}, rbr_pressures{}, rbr_depths_adjusted{}

% For Botany Bay Data
start_time = datenum(2018,6,20,12,0,0);
end_time = datenum(2018,7,26,0,0,0);
cmab = [0.3 0.15]; % RBR 1 & RBR5, which is guessed
load bb_data_RBR1_RBR5.mat
rbr_depths_adjusted = rbr_depths; % Because data not corrected. 
    
% For FIRST Sensor Set
% start_time = datenum(2019,6,4,12,0,0);
% end_time = datenum(2019,7,16,0,0,0);
%        LL  PN  PS  SB  WB  SL   TB
% cmab = [460, 20, 14, 47, 25, 116, 30];

% For SECOND Sensor Set
% load rbr_data_deployment_2.mat
% start_time = datenum(2019,7,19,0,0,0);
% end_time = datenum(2019,8,29,0,0,0);
% % LL, PS, SL, TB, TP, WB
% % N.B. Wall Beach got mangled on Aug 5?try to go before or after this
% % Tomasini only deployed 8/1
% cmab = [460, 14, 116, 30, 3, 25];
% start_time = datenum(2019,8,6,0,0,0);

% Random period to examine
% start_time = datenum(2019,7,26,0,0,0);
% end_time = datenum(2019,7,28,0,0,0);

% Good sea breeze pattern days... 
% start_time = datenum(2019,8,17,0,0,0);
% end_time = datenum(2019,8,24,0,0,0);

% Visible Seiche (LL)...?
% start_time = datenum(2019,8,9,0,0,0);
% end_time = datenum(2019,8,10,0,0,0);

% Peak in Seiching noticeable again (LL?)
% start_time = datenum(2019,7,24,6,0,0);
% end_time = datenum(2019,7,25,18,0,0);

% For THIRD Sensor Set
% load rbr_data_deployment_3.mat
% start_time = datenum(2019,8,30,0,0,0);
% end_time = datenum(2019,9,27,0,0,0);
% % Note big time jump error in LL on Sept 23rd
% % Note Tomasini Point data goes UNTIL 22 NOVEMBER?uses same as 4th deploy
% % LL, PPN, SB, SL, TP
% cmab = [460, 20, 14, 116, 3];
% end_time = datenum(2019,9,22,0,0,0);

% For FOURTH Sensor Set
% start_time = datenum(2019,9,28,0,0,0);
% end_time = datenum(2019,11,21,0,0,0);
% load rbr_data_deployment_4.mat
% % LL, SB, PPN, TP
% % NOTE LAWSONS LANDING TIMESERIES HERE IS MESSED UP, CAME LOOSE SOME WAY IN...
% cmab = [460, 14, 20, 3];
% end_time = datenum(2019,10,20,0,0,0);

% DAYS OF INTEREST
% Sept 28-29 - Big Wind Event, Spring Tide, No particular swell
% start_time = datenum(2019,9,28,12,0,0);
% end_time = datenum(2019,9,29,12,0,0);

% % Oct 19-Oct 19 - Windy, swell arrives, some
% start_time = datenum(2019,10,17,0,0,0);
% end_time = datenum(2019,10,19,12,0,0);
% 
% % Nov 15-16 - Windy, bigger swell arrives, Lukas' Birthday
% start_time = datenum(2019,11,15,6,0,0);
% end_time = datenum(2019,11,16,18,0,0);

% % Nov 20-21 - bigger waves, from local wind? 
% start_time = datenum(2019,11,19,0,0,0);
% end_time = datenum(2019,11,20,12,0,0);

% Nov 7-8 - 12s consistent swell, very light winds
% start_time = datenum(2019,11,7,12,0,0);
% end_time = datenum(2019,11,8,18,0,0);

%% Control Choices

% ONLY USE THE NEXT LINE IF YOU ARE NOT MODIFYING SENSOR_CHOICE OUTSIDE
% THIS SCRIPT (I.E. IN META_RUN.M)
sensor_choice = 1;

% Mechanical Options
% "1" turns them on, "0" turns them off
adjust_pressure = 1; % Default is 1
variance_preserving = 1; % Default is 1
use_t_tide = 1; % Default is 1, turn off if very short time window...
use_windowing = 0; % Default is 0
use_median_power = 0; % Not fully-functional yet. 

% Visualization Options
% "1" turns them on, "0" turns them off
make_spectra_plot = 1;
visualize_conditions = 1;
use_residual_spectra = 0; 
do_energy_wave_correlations = 1;
visualize_big_spectra = 1;

% Spectrum Choice: 0 for none
%                  1 for High/Trans/Low Tide
%                  2 for Morning/Midday/Evening
%                  3 for Flooding/Ebbing/Slack
spectra_choice = 3;

% ALL THREE OF THE NEXT VARIABLES ARE IN HOURS
ea_spacing = 1; % What is the spacing between each ensemble?
window_length = 3; % Each ensemble represents how many hours? 
instance_length = 0.75; % What is the length of each instance in an ensemble?

% 1 /  3 / 0.75 has settled out to be the default
% 1 /  6 / 3 is useful for capturing seiche stuff

%% Setup

% Cutoffs between types of waves (in seconds)
max_period_wind = 4;
max_period_swell = 25;
max_period_igw = 250;

% How many points do you want to use for the moving average of S? 
n_smooth = 5;

% max_varpreserv_power = 1*10^-4; % good for TB data
max_varpreserv_power = 3*10^-3; % good for BB data

fprintf(['*** RUNNING CODE FOR ' labels{sensor_choice} ' with %.2f, %.2f, %.2f. ***\n'],ea_spacing,window_length,instance_length);

if end_time - start_time < 0.75
    use_t_tide = 0;
    fprintf('Time window less than 18 hours: will not use t_tide to de-tide.\n');
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

times = datenum(rbr_times{sensor_choice});
start_index = find(times == start_time);
end_index = find(times == end_time);
if isempty(start_index) || isempty(end_index)
    error('Start or End Time index not found. Check dates and such.');
end
trimmed_depth_signal = depth_signal(start_index:end_index); 
trimmed_times = rbr_times{sensor_choice}(start_index:end_index); % Comes as datetimes

% Remove the tidal signal using t_tide (https://www.eoas.ubc.ca/~rich/)
if use_t_tide
    [~,d_out] = t_tide(trimmed_depth_signal,T/3600,'output','none');
    eta = trimmed_depth_signal - d_out;
else
    eta = trimmed_depth_signal;
    % Unadjusted for tide!!!
end

% Flatten any bottoming-out chunks ...?
if min(trimmed_depth_signal) <= 0, fprintf('!!! Depth signal bottoms out. Take note! Setting values to 0 !!!\n'), end    
trimmed_depth_signal(trimmed_depth_signal <= 0.01) = 0;

mab = cmab./100;
tide_range = max(trimmed_depth_signal) - min(trimmed_depth_signal); % in m

%% Machinery

[W,big_average_S] = f_ensemble_average_spectra(T,eta,instance_length*60*60,5);

if adjust_pressure
    big_average_S = Sd_to_Ss(trimmed_depth_signal,mab(sensor_choice),W,big_average_S);
end

measurements_per_window = window_length*60*60*fs;

% Decision: EA spacing dominates, we will clip timeseries at end if the
% combination of parameters doesn't "fit" perfectly into the timeseries. 

n_windows = floor((length(eta) - measurements_per_window)/(ea_spacing*60*60*fs));

matrixS = zeros(n_windows,length(big_average_S));
matrixSfreq = zeros(size(matrixS));
% Each row in matrix_S is a power spectra for window. 
% The earlier time is at the top, the latter time at the bottom. 

window_depths = zeros(1,n_windows);
window_slope = zeros(1,n_windows);
window_Hs_wind = zeros(1,n_windows);
window_Hs_swell = zeros(1,n_windows);
window_Hs_igw = zeros(1,n_windows);
window_Hs_seiche = zeros(1,n_windows);
window_times = [datetime(2020,1,1,1,1,1)]; % To initiate it
running_S = zeros(size(big_average_S));

high_tide_running_S = zeros(1,length(big_average_S));
n_high_tide_S = 0;
med_tide_running_S = zeros(1,length(big_average_S));
n_med_tide_S = 0;
low_tide_running_S = zeros(1,length(big_average_S));
n_low_tide_S = 0;

morning_running_S = zeros(1,length(big_average_S));
n_morning_S = 0;
afternoon_running_S = zeros(1,length(big_average_S));
n_afternoon_S = 0;
evening_running_S = zeros(1,length(big_average_S));
n_evening_S = 0;

flooding_running_S = zeros(1,length(big_average_S));
n_flooding_S = 0;
ebbing_running_S = zeros(1,length(big_average_S));
n_ebbing_S = 0;
slack_running_S = zeros(1,length(big_average_S));
n_slack_S = 0;

for nn = 0:n_windows-1
    si = nn*ea_spacing*60*60*fs + 1;
    ei = nn*ea_spacing*60*60*fs + measurements_per_window;
    
    if ei > length(eta), error('Encountered an issue with indexing the signal.'), end
    
    ensemble = eta(si:ei);
    ensemble = detrend(ensemble,2); % 2 or 3?
    
    [fitcoeff,~,~] = polyfit(datenum(trimmed_times(si:ei)),trimmed_depth_signal(si:ei),1);
    window_slope(nn+1) = fitcoeff(1);
    
    window_depths(nn+1) = mean(trimmed_depth_signal(si:ei));
    
    if use_windowing
        % Uses Tukey Window, per recommendation here: https://arena-attachments.s3.amazonaws.com/5331399/0509f4edab72d8fb7e4628a7f6fcea3a.pdf?1571926476
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
            if mean(trimmed_depth_signal(si:ei)) > min(trimmed_depth_signal)+(2/3)*tide_range
                high_tide_running_S = high_tide_running_S + S;
                n_high_tide_S = n_high_tide_S + 1;
            elseif mean(trimmed_depth_signal(si:ei)) < min(trimmed_depth_signal)+(1/3)*tide_range
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
    end
    
    matrixS(nn+1,:) = movmean(S,5);
    running_S = running_S + S;
    freq = W./(2*pi);
    matrixSfreq(nn+1,:) = matrixS(nn+1,:).*freq;
    
    % Wind 
    [~,si] = min(abs(freq - 1/max_period_wind));
    m0_wind = trapz(freq(si:end),S(si:end));

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

    window_Hs_wind(nn+1) = 4*sqrt(m0_wind);
    window_Hs_swell(nn+1) = 4*sqrt(m0_swell);
    window_Hs_igw(nn+1) = 4*sqrt(m0_igw);
    window_Hs_seiche(nn+1) = 4*sqrt(m0_seiche);
end

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
end
running_S = running_S./n_windows;


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
m0_wind = trapz(freq(si:end),running_S(si:end));

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

Hs_wind = 4*sqrt(m0_wind);
Hs_swell = 4*sqrt(m0_swell);
Hs_igw = 4*sqrt(m0_igw);
Hs_seiche = 4*sqrt(m0_seiche);

fprintf('Overall Hs wind = %f, Hs swell = %f, Hs IGW = %f, Hs Seiche = %f\n',Hs_wind,Hs_swell,Hs_igw,Hs_seiche);
figure
pie([m0_wind m0_swell m0_igw m0_seiche].*(10^7),{'Wind','Swell','IGW','Seiche'});
title(['Energy (m^2) Budget Partition for ' labels{sensor_choice}]);

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

%% Correlations

if do_energy_wave_correlations

%     % Wind
%     [~,si] = min(abs(freq - 1/max_period_wind));
%     wind_energy = mean(matrixS(:,si:end),2);
%     
%     % Swell
%     [~,si] = min(abs(freq - 1/max_period_swell));
%     [~,ei] = min(abs(freq - 1/max_period_wind));
%     swell_energy = mean(matrixS(:,si:ei),2);
%     
%     % IGW
%     [~,si] = min(abs(freq - 1/max_period_igw)); 
%     [~,ei] = min(abs(freq - 1/max_period_swell));
%     igw_energy = mean(matrixS(:,si:ei),2);
%     
%     % Seiche
%     [~,ei] = min(abs(freq - 1/max_period_igw));
%     seiche_energy = mean(matrixS(:,2:ei),2);

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
%     R = corrcoef(wind_for_corr,wind_energy); R = R(2)
%     corrcoef(sni_for_corr(4:1347),igw_energy(4:1347))

%     [acf,lags,~] = autocorr(igw_energy);
%     plot(lags.*ea_spacing,acf);

    figure
    plot(window_times,window_Hs_wind);
    hold on
    plot(window_times,window_Hs_swell);
    plot(window_times,window_Hs_igw);
    plot(window_times,window_Hs_seiche);
    yyaxis right
    plot(trimmed_times,trimmed_depth_signal,'g--');
    legend('H_s Wind','H_s Swell','H_s IGW','H_s Seiche','Depth');
    title([labels{sensor_choice} ' H_s, ' datestr(start_time) ' to ' datestr(end_time) '. ' extra]);
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
        matrixSfreq(2,end) = max_varpreserv_power;
        matrixSfreq(2,end-1) = 0;
        contourf(datenum(window_times(:)),logfreq(2:end),matrixSfreq(2:end,:),15,'LineColor','none');
    else
        % Impose high & low values to fix scale bar. 
        logSt(2,end) = -1;
        logSt(2,end-1) = -9;
        contourf(datenum(window_times(:)),logfreq(2:end),logSt(2:end,:),15,'LineColor','none'); % 9 is pretty good, or 15
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
            semilogx(freq,median(matrixS).*freq,'k');
        else
            semilogx(freq,running_S.*freq,'k');
        end
        hold on
        
        switch spectra_choice
            case 1
                semilogx(freq,high_tide_running_S.*freq,'Color',colors(1,:));
                semilogx(freq,med_tide_running_S.*freq,'Color',colors(2,:));
                semilogx(freq,low_tide_running_S.*freq,'Color',colors(3,:));
                legend('Total Running','High Tide','Med Tide','Low Tide');
            case 2
                semilogx(freq,morning_running_S.*freq,'r');
                semilogx(freq,afternoon_running_S.*freq,'g');
                semilogx(freq,evening_running_S.*freq,'b');
                legend('Total Running','Morning','Afternoon','Evening');
            case 3
                semilogx(freq,flooding_running_S.*freq,'r');
                semilogx(freq,ebbing_running_S.*freq,'g');
                semilogx(freq,slack_running_S.*freq,'b');
                legend('Total Running','Flooding','Ebbing','Slack');
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
        loglog(freq,running_S,'b');
        axis([min(freq),max(freq),10^-8,10^-2]);
        ylabel('Depth Power Density (m^2/Hz)');
    end
    xlabel('Frequency (Hz)');
    
    % OK actually we're just gonna fit to the juicy middle part...
    [~,si] = min(abs(freq - 10^-2));
    [~,ei] = min(abs(freq - 10^-1));
    
    % Or maybe just 0.05 Hz to the end, per Hughes et al 2014
%     [~,si] = min(abs(freq - 0.05)); % to end

    fitcoeff = polyfit(logfreq(si:ei),logBigS(si:ei),1);
    
    fprintf('Spectra curve slope is %f\n',fitcoeff(1));

    title(['Mean of ' num2str(instance_length) '-hr Spectra at ' labels{sensor_choice} ' from ' datestr(start_time) ' to ' datestr(end_time)]);
end

%% Energy Flux & Dissipation Stuff

% First, need to find energy flux for each frequency bucket. 
% Energy Flux = Group Celerity * Energy Density 
% (really needs to then be multiplied by rho*g*(1/8)...

g = 9.80665;

% All of these are vectors, values for each frequency (really, W) value
kh = qkhfs(W,mean(trimmed_depth_signal)); % Oy, taking in a lot of time here...
k = kh./mean(trimmed_depth_signal); % a bit sloppy
Cp = ((g./k).*tanh(kh)).^0.5;
Cg = Cp.*(0.5+kh./sinh(2.*kh));

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
