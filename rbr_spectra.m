% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last Edited 10 Dec 2019

% load tomales_rbrs.mat
% tomales_rbrs has labels{}, rbr_depths{}, rbr_times{}, rbr_pressures{}, rbr_depths_adjusted{}

n_sensors = length(labels);
fs = 2; % Hz
T = 1/fs; % in seconds

% figure

for kk = 1:n_sensors
    depth_signal = rbr_depths_adjusted{kk};
    times = datenum(rbr_times{kk}); % Converted to datenum for easy working
    site = labels{kk};
    
    % July 10 - Low Winds, long-period swell
    start_index = find(times == datenum(2019,7,10,0,0,0));
    end_index = find(times == datenum(2019,7,11,0,0,0));
    
    % July 4 - High Winds, long-period swell ish
%     start_index = find(times == datenum(2019,7,4,0,0,0));
%     end_index = find(times == datenum(2019,7,5,0,0,0));
    
    % June 14 - High Winds, only wind chop
%     start_index = find(times == datenum(2019,6,14,0,0,0));
%     end_index = find(times == datenum(2019,6,15,0,0,0));
    
    % June 20 - Low Winds, mid-period swell
%     start_index = find(times == datenum(2019,6,20,0,0,0));
%     end_index = find(times == datenum(2019,6,21,0,0,0));
% 
    % June 27 - Low winds, 14s swell
%     start_index = find(times == datenum(2019,6,27,0,0,0));
%     end_index = find(times == datenum(2019,6,28,0,0,0));
    
    snippet = depth_signal(start_index:end_index);
    
    %%%%% REGULAR FFT
%     [W,E] = ezfft(T,snippet);
    
    %%%%% FREQUENCY AVERAGING BEFORE ENSEMBLE AVERAGING
%     E_fa = movmean(E,13); 
     
    %%%%% ENSEMBLE AVERAGING
    [W,E_ea] = f_ensemble_average_spectra(T,snippet,60*60*3); % 3 hours
    
    %%%%% FREQUENCY AVERAGING AFTER ENSEMBLE AVERAGING
    freq = W./(2*pi); 
    E_ea_fa = movmean(E_ea,5);
    
    % Normal Loglog plots
    loglog(freq,E_ea_fa);
    
    % Variance-Preserving Spectra
%     figure
%     semilogx(freq, E_ea.*freq);
    
%     xlabel('Frequency (Hz)');
%     ylabel('Spectra Power');
%     ylabel('Variance per cycle'); % Use if Variance-Preserving Spectra
%     title(labels{kk})
    hold on
end

% Infragravity Range... 30s to 5min (Munk)

legend(labels);
xlabel('Frequency (Hz)');
ylabel('Spectral Power');
% ylabel('Variance per cycle'); % Use if Variance-Preserving Spectra

title('Hourlong Ensemble Averaging, Frequency Averaged on July 10');
% title('Half Hour Ensemble Avg, Low Low on 28 June, 2-5:00');
