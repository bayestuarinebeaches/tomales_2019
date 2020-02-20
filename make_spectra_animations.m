% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last Edited 4 November 2019

get_ocean_conditions
load tomales_rbrs.mat

n_sensors = length(labels);
fs = 2; % Hz
T = 1/fs; % in seconds
order = [1 5 6 3 4 8 7]; % Order of sites from North to South

for kk = 1:n_sensors
    depth_signal = rbr_depths_adjusted{kk};
    times = datenum(rbr_times{kk}); % Converted to datenum for easy working

    %%%%% INPUT YOUR DESIRED TIME WINDOW HERE ...
%     start_time = datenum(2019,7,15,23,0,0);
%     end_time = datenum(2019,7,16,7,0,0);
    start_time = datenum(2019,7,15,20,0,0);
    end_time = datenum(2019,7,16,2,0,0);
    
    n_hours = (end_time-start_time)*24; % Because datenums are in days.
    start_index = find(times == start_time);
    end_index = find(times == end_time);
    
    snippet = depth_signal(start_index:end_index);
%     snippet = detrend(snippet); % Caution, makes the signal not look normal
    
    gif_title = ['test_',labels{kk},'.gif'];
    gif_delay = 0.5; % seconds
    
    % STRONG TO SOFT RED/BLUE COLORING
    start_color_r = [255,0,0]./255;
    end_color_r = [255,200,200]./255;
    colors_r = [linspace(start_color_r(1),end_color_r(1),n_hours+1)',linspace(start_color_r(2),end_color_r(2),n_hours+1)',linspace(start_color_r(3),end_color_r(3),n_hours+1)'];
    start_color_b = [0,0,255]./255;
    end_color_b = [200,200,255]./255;
    colors_b = [linspace(start_color_b(1),end_color_b(1),n_hours+1)',linspace(start_color_b(2),end_color_b(2),n_hours+1)',linspace(start_color_b(3),end_color_b(3),n_hours+1)'];
    
    %%%%% ENSEMBLE AVERAGING
    measurements_per_hour = fs*60*60; 
    running_E = zeros(1,measurements_per_hour/2);
    subplot(4,2,order(kk))
    for hh = 0:(n_hours-1) % Adjust this to go through all intervals in your snippet. 
        [W,E] = ezfft(T,snippet((hh*measurements_per_hour+1):(hh+1)*measurements_per_hour),'hann');
        running_E = running_E + E;
        freq = W./(2*pi); period = 1./freq; % in Seconds
        
%         spot_in_tidal_range = round((snippet(hh*measurements_per_hour+1)-min(snippet))/(max(snippet)-min(snippet))*n_hours)+1;
        
%         figure
%         subplot(3,2,1:2);
%         plot(rbr_times{kk}(start_index:end_index),snippet,'b-');
%         hold on
%         scatter(rbr_times{kk}(start_index+hh*measurements_per_hour+1),snippet(hh*measurements_per_hour+1),'r*');
%         ylabel('Depth (m)');
%         title('Water Depth at Sensor');

        %%%%% FREQUENCY AVERAGING
        upper_cutoff = find(abs(freq - 0.003) < 0.0001);
        igw_cutoff = find(abs(freq-0.010) < 0.0001); % Using 0.05Hz as IGW cutoff
        E = [movmean(E(1:upper_cutoff),3) movmean(E(upper_cutoff+1:igw_cutoff),7) movmean(E(igw_cutoff+1:end),17)];
        
        % Normal Loglog plots
%         subplot(3,2,3:6);
%         loglog(freq,E);
%         ylim([10^-11 1.1]);
%         xlim([10^-4 10^0]);
%         xlabel('Frequency'); ylabel('Spectral Power');
% 
%         title(['Frequency-Averaged Spectra ',datestr(start_time),' ',num2str(hh),'th hour at ',labels{kk}]);
%         saveas(gcf,'temp.png');
%         [A,~] = imread('temp.png');
%         [X,map] = rgb2ind(A,256);
%         if hh == 0
%             imwrite(X,map,gif_title,'gif','LoopCount',inf,'DelayTime',gif_delay);
%         else
%             imwrite(X,map,gif_title,'gif','WriteMode','append','DelayTime',gif_delay);
%         end

%         loglog(freq,E,'Color',colors_b(hh+1,:));
%         semilogx(freq,E.*freq,'Color',colors_r(spot_in_tidal_range,:)); %
%         hold on
        
    end
    
    E = running_E./n_hours; % Divide by 24 hours in a day to finish EA
% %     
% %     %%%%% FREQUENCY AVERAGING
      upper_cutoff = find(abs(freq - 0.003) < 0.0001);
      igw_cutoff = find(abs(freq-0.010) < 0.0001); % Using 10s as cutoff, so we can keep looking at swell...
      E = [movmean(E(1:upper_cutoff),3) movmean(E(upper_cutoff+1:igw_cutoff),5) movmean(E(igw_cutoff+1:end),17)];
% % 
% %     %%%%% VARIANCE-PRESERVING SPECTRA
% % %     semilogx(freq, E_ea.*freq);
% % %     ylabel('Variance per Spectra);
% %     
% %     %%%%% NORMAL SPECTRA
    loglog(freq,E);

    ylim([10^-11 1.1]);
    xlim([10^-4 10^0]);
    xlabel('Frequency (Hz)'); ylabel('Spectral Power');
    title(labels{kk});
% close all
end

% sgtitle(['Frequency-Averaged Spectra ',datestr(start_time),' ',num2str(n_hours),' hours HHW -> LLW Tide']);
sgtitle(['EA/FA Spectra ',datestr(start_time),' - ',num2str(n_hours),' Hrs Over High Tide']);
% legend('High Tide','','','','','','Low Tide');
