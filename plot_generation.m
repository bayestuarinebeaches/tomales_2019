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
            scatter(tbb_wind.time,tbb_wind.spd.*0.44704,'.');
        case 3
            scatter(hioc_wind.time,hioc_wind.spd,'.');
        case 4
            scatter(tboc_wind.time,tboc_wind.spd,'.');
        case 5
            scatter(tbb_wind.time,tbb_wind.spd.*0.44704,'.');
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

    xlim([window_times(1) window_times(end)]);
    linkaxes(ff,'x');
    xlim([window_times(1) window_times(end)]);
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
        contourf(datenum(window_times(:)),logfreq(n_leakage_ignore:end),log10(matrixSfreq(n_leakage_ignore:end,:)),15,'LineColor','none');
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
        c.Label.String = 'Log_{10} Energy Density (m^2/Hz)';
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
            scatter(datenum(tbb_wind.time),tbb_wind.spd.*0.44704,'.');
        case 3
            scatter(datenum(hioc_wind.time),hioc_wind.spd,'.');
        case 4
            scatter(datenum(tboc_wind.time),tboc_wind.spd,'.');
        case 5
            scatter(tbb_wind.time,tbb_wind.spd.*0.44704,'.');
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
            scatter(bml_wind.dir_time,bml_wind.dir,'+');
    end   
    scatter(swell.time,swell.dir,'go');
    ylabel('Direction (°)');
    ylim([0 360]); % weird outliers sometimes...
    xlim([window_times(1) window_times(end)]);

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
    xlim([window_times(1) window_times(end)]);

%     linkaxes(ff,'x');
    
end

if make_conditions_plot
    figure
    sgtitle('Conditions During Deployment');
    
    ff(1) = subplot(3,2,[1,2]);
    yyaxis left
    ylabel('Wind Speed (m/s)');
    switch wind_option
        case 1
            scatter(bml_wind.spd_time,bml_wind.spd,'.');
        case 2
            scatter(tbb_wind.time,tbb_wind.spd.*0.44704,'.');
        case 3
            scatter(hioc_wind.time,hioc_wind.spd,'.');
        case 4
            scatter(tboc_wind.time,tboc_wind.spd,'.');
        case 5
            scatter(tbb_wind.time,tbb_wind.spd.*0.44704,'.');
    end
    
    ylim([0 10]);
    ylabel('Wind Speed (m/s)');
    yyaxis right
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
    ylabel('Wind Direction (°)');
    ylim([0 360]); % weird outliers sometimes...
    
    ff(2) = subplot(3,2,[3,4]);
    yyaxis left
    scatter(swell.time,swell.hgt,'.');
    ylabel('Offshore Significant Wave Height (m)');
    
    yyaxis right
    hold on
    scatter(swell.time,swell.dir,'+');
    ylabel('Offshore Dominant Wave Direction (°)');
    
    ff(3) = subplot(3,2,[5,6]);
    yyaxis left
    plot(trimmed_times,trimmed_depth_signal,'-');
    ylabel(['Depth from Sensor at ' labels{sensor_choice} ' (m)']);
    
    yyaxis right
    hold on
    scatter(swell.time,swell.per,'.');
    ylabel('Offshore Dominant Wave Period (s)');
    
    xlim([window_times(1) window_times(end)]);
    linkaxes(ff,'x');
    xlim([window_times(1) window_times(end)]);
    
end

wave_direction_mapped = interp1(datenum(swell.time),double(swell.dir),datenum(window_times'));
wave_period_mapped = interp1(datenum(swell.time),swell.per,datenum(window_times'));
wave_height_mapped = interp1(datenum(swell.time),swell.hgt,datenum(window_times'));

switch wind_option
    case 1
        wind_speed_mapped = interp1(datenum(bml_wind.spd_time),bml_wind.spd,datenum(window_times'));
        wind_dir_mapped = interp1(datenum(bml_wind.dir_time),double(bml_wind.dir),datenum(window_times'));
    case 2
        wind_speed_mapped = interp1(datenum(tbb_wind.time),tbb_wind.spd,datenum(window_times')).*0.44704;
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
        wind_speed_mapped = interp1(datenum(tbb_wind.time),tbb_wind.spd,datenum(window_times')).*0.44704;
        wind_dir_mapped = interp1(datenum(bml_wind.dir_time(unique_index)),bml_wind.dir(unique_index),datenum(window_times'));
end
wind_direction_NW = wind_dir_mapped > 280; % in degrees, 4.89 for radians
eval(['Hs_wind_timeseries = Hs_wind_' num2str(sensor_choice) ';']);
eval(['Hs_swell_timeseries = Hs_swell_' num2str(sensor_choice) ';']);
eval(['Hs_igw_timeseries = Hs_igw_' num2str(sensor_choice) ';']);
Hs_wind_timeseries = Hs_wind_timeseries';
Hs_swell_timeseries = Hs_swell_timeseries';
Hs_igw_timeseries = Hs_igw_timeseries';


%% Polar Plotting
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

%% Etc Plotting
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
            scatter(tbb_wind.time,tbb_wind.spd.*0.44704,'.');
        case 3
            scatter(hioc_wind.time, hioc_wind.spd,'.');
        case 4
            scatter(tboc_wind.time, tboc_wind.spd,'.');
        case 5
            scatter(tbb_wind.time,tbb_wind.spd.*0.44704,'.');
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
            scatter(datenum(tbb_wind.time),tbb_wind.spd.*0.44704,'.');
        case 3
            scatter(datenum(hioc_wind.time),hioc_wind.spd,'.');
        case 4
            scatter(datenum(tboc_wind.time),tboc_wind.spd,'.');
        case 5
            scatter(datenum(tbb_wind.time),tbb_wind.spd.*0.44704,'.');
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