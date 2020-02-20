% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last edited 23 October 2019

% Round One
% prepend = 'RBR_data/';
% files = {'20190717_LL_124165/data.txt','20190717_PN_124162/data.txt','20190717_PS_124163/data.txt','20190717_SB_124168/data.txt','20190717_WB_124125/data.txt','20190718_SL_124164/data.txt','20190801_TB_124127/data.txt'};
% labels = {'Lawsons Landing','Pelican Pt N','Pelican Pt S','Seal Beach','Wall Beach','Sacramento Landing','Tomales Beach'};

% Round Two
% prepend = 'RBR_data/20190829/';
% files = {'20190829_LL_124166_data.txt','20190829_PS_124162_data.txt','20190829_SL_124168_data.txt','20190829_TB_124165_data.txt','20190829_TP_124164_data.txt','20190829_WB_124161_data.txt'};
% labels = {'Lawsons Landing','Pelican Pt S','Sacramento Landing','Tomales Beach','Tomasini Point','Wall Beach'};

% Round Three

% Round Four
prepend = 'RBR_data/20191124/';
files = {'20191124_LL_124164_data.txt','20191124_SB_124165_data.txt','20191210_PN_124166_data.txt'};
labels = {'Lawsons Landing','Seal Beach','Pelican Point N'};

fs = 2; % 1/seconds
T = 1/fs;

fprintf('Starting to read RBR Data Files.\n');
rbr_times = {};
rbr_depths = {};
rbr_pressures = {};
rbr_depths_adjusted = {};
for ii = 1:length(files)

    fid = fopen(strcat(prepend,files{ii}),'rt');
    tmp = textscan(fid,'%d-%d-%d %d:%d:%f,%f,%f,%f','Headerlines',1);
    year = tmp{1};
    month = tmp{2};
    day = tmp{3};
    hour = tmp{4};
    minute = tmp{5};
    second = tmp{6};
    pressure = tmp{7};
    % sea_pressure = tmp{8};
    depth = tmp{9};
    
    rbr_time = datetime(tmp{1},tmp{2},tmp{3},tmp{4},tmp{5},tmp{6});
%     density_start_index = find(datenum(density_time) == datenum(rbr_time(1)));
%     pressure_start_index = find(datenum(pressure_time) == datenum(rbr_time(1)));
%     
%     % Use the preceding hourly water density and barometric pressure values to adjust the pressure reading. 
%     for jj = 1:length(pressure)
%         % woof
%         prior_hour = datenum(datetime(tmp{1}(jj),tmp{2}(jj),tmp{3}(jj),tmp{4}(jj),0,0));
%         density_index = find(prior_hour == datenum(density_time));
%         pressure_index = find(prior_hour == datenum(pressure_time));
%     end
    
%     depth = pressure./(1000*9.81);
    
    rbr_times{end+1} = rbr_time;
    rbr_depths{end+1} = depth;
    rbr_pressures{end+1} = pressure;
    
    fprintf('Done with file %d of %d.\n',ii,length(files));
end

get_ocean_conditions

% for kk = 1:length(files)
%     progress_bar(kk,1,length(files));
%     fixed_depths = rbr_pressures{kk};
%     for ll = 1:length(fixed_depths)
%         [~,air_pressure_index] = min(abs(rbr_times{kk}(ll) - pressure_time));
%         [~,water_density_index] = min(abs(rbr_times{kk}(ll) - density_time));
%         % Assumes BOON saltwater density and air pressure represent the
%         % whole bay, which is definitely not true but ok...
%         
%         
%         % rbr_depths are untreated, output from RBR using some assumptions. 
%         % rbr_pressures are direct pressure readings from the sensor. 
%         % Here, we'll adjust pressure values to find fixed depths. 
%         
%         % This line takes the old pressure value (in decibars), converts to
%         % Pascals, subtracts the barotropic pressure, then divides by
%         % density and gravitational acceleration to find height of the
%         % water column above the sensor in meters. 
%         
%         fixed_depths(ll) = (fixed_depths(ll)*10000-bpressure(air_pressure_index)*100)./(water_density(water_density_index)*9.81);
%     end
%     rbr_depths_adjusted{end+1} = fixed_depths;
% end

save tomales_rbrs.mat labels rbr_times rbr_depths rbr_pressures %rbr_depths_adjusted