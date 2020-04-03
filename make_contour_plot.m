% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last Edited 3 April 2020

%% Binary Choices
% "1" turns them on, "0" turns them off
attenuate_pressure = 1;
visualize_conditions = 1;
use_residual_spectra = 0;
use_windowing = 1;
detrend_ensemble = 1; % use esp. if instance length is shorter than tidal period
do_energy_correlations = 1;
visualize_big_spectra = 1;
variance_preserving = 0;

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

sensor_choice = 3;

%% Setup

if instance_length > window_length, error('Ensemble Instance length must be smaller than the averaging window length.'), end
if window_length < ea_spacing, error('Right now, code requires windows be longer than ensemble spacing.'), end
if attenuate_pressure
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

if isempty(start_index), error('Start time not found.'), end

mab = cmab./100;
    

%% Machinery

[W,big_average_S] = f_ensemble_average_spectra(T,trimmed_depth_signal,instance_length*60*60,5);

if attenuate_pressure
    big_average_S = Sp_to_Ss(trimmed_depth_signal,mab(sensor_choice),W,big_average_S);
end

measurements_per_window = window_length*60*60*fs;

% Decision: EA spacing dominates, we will clip timeseries at end if the
% combination of parameters doesn't "fit" perfectly into the timeseries. 

n_insts = floor((length(trimmed_depth_signal) - measurements_per_window)/(ea_spacing*60*60*fs));

matrixS = zeros(n_insts,length(big_average_S));
% Each row in matrix_E is a power spectra for window. 
% The earlier time is at the top, the latter time at the bottom. 

ea_times = [datetime(2020,1,1,1,1,1)]; % To initiate it
runningS = zeros(size(big_average_S));

for nn = 0:n_insts-1
    si = nn*ea_spacing*60*60*fs + 1;
    ei = nn*ea_spacing*60*60*fs + measurements_per_window;
    
    if ei > length(trimmed_depth_signal), error('Encountered an issue with indexing the signal.'), end
    
    if detrend_ensemble
        ensemble = detrend(trimmed_depth_signal(si:ei),1);
        % Detrend type 1 removes a linear trend. 
    else
        ensemble = trimmed_depth_signal(si:ei);
    end
    
    if use_windowing
        % Uses Tukey Window, per recommendation here: https://arena-attachments.s3.amazonaws.com/5331399/0509f4edab72d8fb7e4628a7f6fcea3a.pdf?1571926476
        [W,S] = f_ensemble_average_spectra(T,ensemble.*tukeywin(length(ensemble)),instance_length*60*60,n_smooth);
    else
        [W,S] = f_ensemble_average_spectra(T,ensemble,instance_length*60*60,n_smooth);
    end
     
    ea_times(nn+1) = trimmed_times(nn*ea_spacing*60*60*fs+1+floor(measurements_per_window/2)); % These times are at center of interval.
    
    if attenuate_pressure
        S = Sd_to_Ss(ensemble,mab(sensor_choice),W,S);
    end
    
    matrixS(nn+1,:) = movmean(S,5);
    runningS = runningS + S;
    freq = W./(2*pi);
end

runningS = runningS./n_insts;


%% Adjustment 
    
logfreq = log10(freq);
logS = log10(matrixS);
logSt = zeros(size(logS));
logBigS = log10(big_average_S);
logRunningS = log10(runningS);
if use_residual_spectra
    for nn = 1:n_insts
        logSt(nn,:) = logS(nn,:) - logBigS; % Could subtract logRunning_S instead
    end
else
    logSt = logS;
end

%% Conditions

% get_conditions
load conditions.mat
load tbb_wind.mat

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



%% Correlations

if do_energy_correlations

    % Let local wind waves be 1s to 8s (1 to .125 frequency)
    
    [~,si] = min(abs(freq - 0.125));
    wind_energy = mean(matrixS(:,si:end),2);
    
    % Let swell be 8s to 20s (0.125 to 0.05 frequency)
    
    [~,si] = min(abs(freq - 0.05));
    [~,ei] = min(abs(freq - 0.125));
    swell_energy = mean(matrixS(:,si:ei),2);
    
    % Let igw be 20s to 5min (0.05 to 0.004 0.0033333... frequency)
    
    [~,si] = min(abs(freq - 0.004)); %0.00333));
    [~,ei] = min(abs(freq - 0.05));
    igw_energy = mean(matrixS(:,si:ei),2);
    
    % Let seiches be 5min to 30min (0.003333... to 0.00055555... frequency)
    
    [~,si] = min(abs(freq - 0.000555));
    [~,ei] = min(abs(freq - 0.00333));
    seiche_energy = mean(matrixS(:,si:ei),2);

    % Interpolate conditions onto energy times
    % Note, these don't account for direction of wind or swell. 
    wind_for_corr = interp1(datenum(tbb_wind.time),tbb_wind.spd,datenum(ea_times)); % Note, using TBB wind here
    swell_for_corr = interp1(datenum(swell.time),swell.hgt,datenum(ea_times));
    sni_for_corr = interp1(datenum(sni.time),sni.value,datenum(ea_times));
    
    figure
    hold on
    plot(datenum(ea_times),wind_energy./max(wind_energy));
    plot(datenum(ea_times),swell_energy./max(swell_energy));
    plot(datenum(ea_times),igw_energy./max(igw_energy));
%     plot(datenum(ea_times),seiche_energy./max(seiche_energy));
    plot(datenum(trimmed_times),trimmed_depth_signal./max(trimmed_depth_signal),'g:');
    xlim([datenum(ea_times(1)) datenum(ea_times(end))]);
    title(['Normalized Energy from Wind, Swell, IGW, Seiche Bands at ' labels{sensor_choice}]);
    legend('< 8s','8s to 20s','20s to 250s');%,'5min to 30min');
    
    % For cross- and auto-correlation...
    
%     R = corrcoef(swell_for_corr,swell_energy); R = R(2);
%     corrcoef(sni_for_corr(4:1347),igw_energy(4:1347))

%     [acf,lags,~] = autocorr(igw_energy);
%     plot(lags.*ea_spacing,acf);

end

%% Plotting

logSt = logSt'; % Transpose! 

figure
sgtitle([labels{sensor_choice} ' Spectra, ' num2str(ea_spacing) '-hr Spacing, ' num2str(window_length) '-hr Window, ' num2str(instance_length) '-hour Ensembles ' extra]);

if visualize_conditions
    ff(1) = subplot(4,2,[1,2,3,4]);
end

% It is good to start at 2 to avoid inf. 
contourf(datenum(ea_times(:)),logfreq(2:end),logSt(2:end,:),15,'LineColor','none'); % 9 is pretty good, or 15
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

if visualize_big_spectra
    figure
    if variance_preserving
        semilogx(freq,big_average_S.*freq,'b');
        ylabel('Variance Per Hz');
    else
        loglog(freq,big_average_S,'b');
        ylabel('Depth Power Density (m^2/Hz)');
    end
    xlabel('Frequency (Hz)');

    fitcoeff = polyfit(logfreq(2:end),logBigS(2:end),1);
    fprintf('Spectra curve slope is %f\n',fitcoeff(1));

    title(['Arithmetic Mean of ' num2str(instance_length) '-hour Instance Spectra at ' labels{sensor_choice} ' from ' datestr(start_time) ' to ' datestr(end_time)]);
end


%% To-Do

% - Automate cmab stuff
% - Quiver plot for winds & swell ? 

