% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last Edited 26 April 2020

%% Binary Choices
% "1" turns them on, "0" turns them off
adjust_pressure = 0;
make_spectra_plot = 1;
visualize_conditions = 1;
use_residual_spectra = 0;
use_windowing = 0;
do_energy_correlations = 0;
visualize_big_spectra = 1;
variance_preserving = 1;

%% Data Controls

% ALL THREE OF THE NEXT VARIABLES ARE IN HOURS
ea_spacing = 1; % What is the spacing between each ensemble?
window_length = 3; % Each ensemble represents how many hours? 
instance_length = 0.75; % What is the length of each instance in an ensemble?

% 6/12/4 is a good combo (smear over a tide)
% 3/24/3 can be nice for slow stuff

% 1 /  3 / 0.75
% 3 / 25 / 0.75

% How many points do you want to use for the moving average of S? 
n_smooth = 5;
    
% For FIRST Sensor Set
% start_time = datenum(2019,6,4,12,0,0);
% end_time = datenum(2019,7,16,0,0,0);
% load RBR_data/20190717/tomales_rbrs.mat
% LL, PPN, PPS, SB, WB, SL, TB
cmab = [460, 20, 14, 47, 25, 116, 30];

% For SECOND Sensor Set
% start_time = datenum(2019,7,18,0,0,0);
% end_time = datenum(2019,8,29,0,0,0);
% load RBR_data/20190829/tomales_rbrs.mat
% LL, PPS, SL, TB, TP, WB
% cmab = [460, 14, 116, 30, 3, 25];

% For THIRD Sensor Set
% start_time = datenum(2019,8,30,0,0,0);
% end_time = datenum(2019,9,26,0,0,0);
% load RBR_data/20190927/tomales_rbrs.mat
% NOTE LAWSONS LANDING TIMESERIES HERE IS MESSED UP ...
% LL, PPN, SB, SL
% cmab = [460, 20, 14, 116];

% For FOURTH Sensor Set
start_time = datenum(2019,9,28,0,0,0);
end_time = datenum(2019,11,24,0,0,0);
% load RBR_data/20191124/tomales_rbrs.mat
load rbr_data_deployment_4.mat
% LL, SB, PPN
% Note LL data in this deployment is questionable
cmab = [460, 14, 20];

% tomales_rbrs.mat has labels{}, rbr_depths{}, rbr_times{}, rbr_pressures{}, rbr_depths_adjusted{}

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
% 
% % Nov 20-21 - bigger waves, from local wind? 
% start_time = datenum(2019,11,19,0,0,0);
% end_time = datenum(2019,11,20,12,0,0);

% start_time = datenum(2019,11,19,12,0,0);
% end_time = datenum(2019,11,20,12,0,0);

sensor_choice = 2;

%% Setup

fprintf(['Running stuff for ' labels{sensor_choice} 'with %d, %d, %d.\n'],ea_spacing,window_length,instance_length);

if instance_length > window_length, error('Ensemble Instance length must be smaller than the averaging window length.'), end
if window_length < ea_spacing, error('Right now, code requires windows be longer than ensemble spacing.'), end
if adjust_pressure
    extra = 'with attenuation';
else
    extra = '';
end

T = seconds(rbr_times{1}(2)-rbr_times{1}(1)); % seconds
fs = 1/T; % Hz

depth_signal = rbr_depths_adjusted{sensor_choice};

times = datenum(rbr_times{sensor_choice});
start_index = find(times == start_time);
end_index = find(times == end_time);
trimmed_depth_signal = depth_signal(start_index:end_index); 
trimmed_times = rbr_times{sensor_choice}(start_index:end_index); % Comes as datetimes

% Remove the tidal signal using t_tide (https://www.eoas.ubc.ca/~rich/)
[tidestruc,d_out] = t_tide(trimmed_depth_signal,T/3600,'output','none');
eta = trimmed_depth_signal - d_out;

if isempty(start_index), error('Start time not found.'), end

mab = cmab./100;
tide_range = max(trimmed_depth_signal) - min(trimmed_depth_signal);

% Cutoffs between types of waves (in seconds)
max_period_wind = 2.5; 
max_period_swell = 25;
max_period_igw = 250;

%% Machinery

[W,big_average_S] = f_ensemble_average_spectra(T,eta,instance_length*60*60,5);

if adjust_pressure
    big_average_S = Sd_to_Ss(eta,mab(sensor_choice),W,big_average_S);
end

measurements_per_window = window_length*60*60*fs;

% Decision: EA spacing dominates, we will clip timeseries at end if the
% combination of parameters doesn't "fit" perfectly into the timeseries. 

n_insts = floor((length(eta) - measurements_per_window)/(ea_spacing*60*60*fs));

matrixS = zeros(n_insts,length(big_average_S));
matrixSfreq = zeros(size(matrixS));
% Each row in matrix_E is a power spectra for window. 
% The earlier time is at the top, the latter time at the bottom. 
high_tide_running_S = zeros(1,length(big_average_S));
n_high_tide_S = 0;
med_tide_running_S = zeros(1,length(big_average_S));
n_med_tide_S = 0;
low_tide_running_S = zeros(1,length(big_average_S));
n_low_tide_S = 0;

ea_times = [datetime(2020,1,1,1,1,1)]; % To initiate it
running_S = zeros(size(big_average_S));

for nn = 0:n_insts-1
    si = nn*ea_spacing*60*60*fs + 1;
    ei = nn*ea_spacing*60*60*fs + measurements_per_window;
    
    if ei > length(eta), error('Encountered an issue with indexing the signal.'), end
    
    ensemble = eta(si:ei);
    ensemble = detrend(ensemble,3);
    
    if use_windowing
        % Uses Tukey Window, per recommendation here: https://arena-attachments.s3.amazonaws.com/5331399/0509f4edab72d8fb7e4628a7f6fcea3a.pdf?1571926476
        [W,S] = f_ensemble_average_spectra(T,ensemble.*tukeywin(length(ensemble)),instance_length*60*60,n_smooth);
    else
        [W,S] = f_ensemble_average_spectra(T,ensemble,instance_length*60*60,n_smooth);
    end
     
    ea_times(nn+1) = trimmed_times(nn*ea_spacing*60*60*fs+1+floor(measurements_per_window/2)); % These times are at center of interval.
    
    if adjust_pressure
        S = Sd_to_Ss(ensemble,mab(sensor_choice),W,S);
    end
    
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
        
    
    matrixS(nn+1,:) = movmean(S,5);
    running_S = running_S + S;
    freq = W./(2*pi);
    
    matrixSfreq(nn+1,:) = matrixS(nn+1,:).*freq;
end

high_tide_running_S = high_tide_running_S./n_high_tide_S;
med_tide_running_S = med_tide_running_S./n_med_tide_S;
low_tide_running_S = low_tide_running_S./n_low_tide_S;
running_S = running_S./n_insts;


%% Adjustment 
    
logfreq = log10(freq);
logS = log10(matrixS);
% logSfreq = log10(matrixSfreq);
logSt = zeros(size(logS));
logBigS = log10(big_average_S);
logRunningS = log10(running_S);
if use_residual_spectra
    for nn = 1:n_insts
        logSt(nn,:) = logS(nn,:) - logRunningS; % Could subtract logBig_S instead, but it doesn't make as much sense with windowing and such. 
    end
else
    logSt = logS;
end

%% Conditions

% get_conditions
load conditions.mat
load tbb_wind.mat
load buoy_spectra.mat

% Spring-Neap Index
% Every m hours, compute index based off min-max...?
% Kinda janky to smooth it out enough but let's give it a shot. 
m = 8;
n_chunks = floor(length(trimmed_depth_signal)./(m*60*60*fs));
diffs = zeros(1,n_chunks);
sni.time = [datetime(2020,1,1,1,1,1)];
for ii = 0:n_chunks-1
    si = floor(ii*m*60*60*fs + 1);
    ei = floor((ii+1)*m*60*60*fs);
    diffs(ii+1) = max(trimmed_depth_signal(si:ei)) - min(trimmed_depth_signal(si:ei));
    sni.time(ii+1) = trimmed_times(ei);
end
tmp_diffs = movmean(movmean(diffs,3),8);
sni.value = (tmp_diffs-min(tmp_diffs))./(max(tmp_diffs)-min(tmp_diffs));

%% Moments

% Wind 
[~,si] = min(abs(freq - 1/max_period_wind));
m0_wind = sum(running_S(si:end).*freq(si:end));

% Swell
[~,si] = min(abs(freq - 1/max_period_swell));
[~,ei] = min(abs(freq - 1/max_period_wind));
m0_swell = sum(running_S(si:ei).*freq(si:ei));

% IGW (0.05 to 0.004 Hz)
[~,si] = min(abs(freq - 1/max_period_igw)); 
[~,ei] = min(abs(freq - 1/max_period_swell));
m0_igw = sum(running_S(si:ei).*freq(si:ei));

Hs_wind = 4*sqrt(m0_wind);
Hs_swell = 4*sqrt(m0_swell);
Hs_igw = 4*sqrt(m0_igw);

fprintf('Hs wind = %f, Hs swell = %f, Hs IGW = %f\n',Hs_wind,Hs_swell,Hs_igw);

%% Correlations

if do_energy_correlations

    % Wind
    [~,si] = min(abs(freq - 1/max_period_wind));
    wind_energy = mean(matrixS(:,si:end),2);
    
    % Swell
    [~,si] = min(abs(freq - 1/max_period_swell));
    [~,ei] = min(abs(freq - 1/max_period_wind));
    swell_energy = mean(matrixS(:,si:ei),2);
    
    % IGW
    [~,si] = min(abs(freq - 1/max_period_igw)); 
    [~,ei] = min(abs(freq - 1/max_period_swell));
    igw_energy = mean(matrixS(:,si:ei),2);
    
    % Seiche
    [~,ei] = min(abs(freq - 1/max_period_igw));
    seiche_energy = mean(matrixS(:,2:ei),2);

    % Interpolate conditions onto energy times
    % Note, these don't account for direction of wind or swell. 
    wind_for_corr = interp1(datenum(tbb_wind.time),tbb_wind.spd,datenum(ea_times)); % Note, using TBB wind here
    wind_for_corr = interp1(datenum(wind.time),wind.spd,datenum(ea_times));
    swell_for_corr = interp1(datenum(swell.time),swell.hgt,datenum(ea_times));
    sni_for_corr = interp1(datenum(sni.time),sni.value,datenum(ea_times));
    
    figure
    hold on
    plot(datenum(ea_times),wind_energy./max(wind_energy));
    plot(datenum(ea_times),swell_energy./max(swell_energy));
    plot(datenum(ea_times),igw_energy./max(igw_energy));
    plot(datenum(ea_times),seiche_energy./max(seiche_energy));
    plot(datenum(trimmed_times),eta./max(eta),'g:');
    xlim([datenum(ea_times(1)) datenum(ea_times(end))]);
    title(['Normalized Energy from Wind, Swell, IGW, Seiche Bands at ' labels{sensor_choice}]);
    legend('Sea','Swell','IGW','Seiche','water surface');
    
    % For cross- and auto-correlation...
    R = corrcoef(wind_for_corr,wind_energy); R = R(2)
%     corrcoef(sni_for_corr(4:1347),igw_energy(4:1347))

%     [acf,lags,~] = autocorr(igw_energy);
%     plot(lags.*ea_spacing,acf);

end

%% Plotting

logSt = logSt'; % Transpose! 
matrixSfreq = matrixSfreq';

if make_spectra_plot
    figure
    sgtitle([labels{sensor_choice} ' Spectra, ' num2str(ea_spacing) '-hr Spacing, ' num2str(window_length) '-hr Window, ' num2str(instance_length) '-hour Ensembles ' extra]);

    if visualize_conditions
        ff(1) = subplot(4,2,[1,2,3,4]);
    end

    % It is good to start at 2 to avoid inf. 
    if variance_preserving
        % Impose high & low values to fix scale bar. 
        matrixSfreq(2,end) = 6*10^-4;
        matrixSfreq(2,end-1) = 0;
        contourf(datenum(ea_times(:)),logfreq(2:end),matrixSfreq(2:end,:),15,'LineColor','none');
    else
        % Impose high & low values to fix scale bar. 
        logSt(2,end) = -1;
        logSt(2,end-1) = -9;
        contourf(datenum(ea_times(:)),logfreq(2:end),logSt(2:end,:),15,'LineColor','none'); % 9 is pretty good, or 15
    end
    c = colorbar('east');
    if use_residual_spectra
        c.Label.String = 'Residual Power Density (m^2/Hz)';
    else
        c.Label.String = 'Power Density (m^2/Hz)';
    end

    % datetick('x',21);
    xlim([datenum(ea_times(1)) datenum(ea_times(end))]);
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

        scatter(datenum(wind.time),wind.spd,'.');
        ylim([0 20]);
        ylabel('Wind Speed (m/s)');
        yyaxis left
        hold on
        scatter(datenum(wind.time),wind.dir,'o');
        scatter(datenum(swell.time),swell.dir,'g+');
        ylabel('Direction (°)');
        ylim([0 360]); % weird outliers sometimes...
    %     xlim([datetime(start_time,'ConvertFrom','datenum') datetime(end_time,'ConvertFrom','datenum')]);
        xlim([start_time end_time]);

        ff(3) = subplot(4,2,[7,8]);
        yyaxis left
        plot(datenum(swell.time),swell.hgt,'+');
        ylabel('H_s (m)');
        hold on
    %     plot(datenum(rbr_times{sensor_choice}),depth_signal,'g-');
        plot(datenum(sni.time),sni.value.*5,'g-');
        plot(datenum(trimmed_times),trimmed_depth_signal,'k-');
        yyaxis right
        scatter(datenum(swell.time),swell.per,'.');
        ylabel('Period (s)');
        ylim([0 20]);
    %     xlim([datetime(start_time,'ConvertFrom','datenum') datetime(end_time,'ConvertFrom','datenum')]);
        xlim([start_time end_time]);

    end

    linkaxes(ff,'x');
end

%% "Big Spectrum" Plotting

if visualize_big_spectra
    figure
    if variance_preserving
        semilogx(freq,running_S.*freq,'k');
        hold on
        semilogx(freq,high_tide_running_S.*freq,'b');
        semilogx(freq,med_tide_running_S.*freq,'r');
        semilogx(freq,low_tide_running_S.*freq,'g');
        ylabel('Variance Per Hz');
        legend('Total Running','High Tide','Med Tide','Low Tide');
        axis([min(freq),max(freq),0,3*10^-5]);
    else
%         loglog(freq,big_average_S,'b');
%         hold on
        loglog(freq,running_S,'b');
        axis([min(freq),max(freq),10^-8,10^-2]);
        ylabel('Depth Power Density (m^2/Hz)');
    end
    xlabel('Frequency (Hz)');
    
    % OK actually we're just gonna fit to the juicy middle part...
    [~,si] = min(abs(freq - 0.001));
    [~,ei] = min(abs(freq - 0.1));
    
    % Or maybe just 0.05 Hz to the end, per Hughes et al 2014
%     [~,si] = min(abs(freq - 0.05));
    
    % Just get the high-freq end?
%     [~,si] = min(abs(freq - 0.4)); % Use this to test the high-freq end

%     fitcoeff = polyfit(logfreq(si:ei),logBigS(si:ei),1);
    fitcoeff = polyfit(logfreq(si:end),logBigS(si:end),1);
    
    fprintf('Spectra curve slope is %f\n',fitcoeff(1));

    title(['Arithmetic Mean of ' num2str(instance_length) '-hour Instance Spectra at ' labels{sensor_choice} ' from ' datestr(start_time) ' to ' datestr(end_time)]);
end


