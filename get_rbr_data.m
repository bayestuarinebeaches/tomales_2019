% Lukas WinklerPrins lukas_wp@berkeley.edu
% Last edited 16 March 2020

% Round One
% prepend = 'RBR_data/20190717/';
% files = {'20190717_LL_124165_data.txt','20190717_PN_124162_data.txt','20190717_PS_124163_data.txt','20190717_SB_124168_data.txt','20190717_WB_124125_data.txt','20190718_SL_124164_data.txt','20190801_TB_124127_data.txt'};
% labels = {'Lawsons Landing','Pelican Pt N','Pelican Pt S','Seal Beach','Wall Beach','Sacramento Landing','Tomales Beach'};

% Round Two
% prepend = 'RBR_data/20190829/';
% files = {'20190829_LL_124166_data.txt','20190829_PS_124162_data.txt','20190829_SL_124168_data.txt','20190829_TB_124165_data.txt','20190829_TP_124164_data.txt','20190829_WB_124161_data.txt'};
% labels = {'Lawsons Landing','Pelican Pt S','Sacramento Landing','Tomales Beach','Tomasini Point','Wall Beach'};

% Round Three
prepend = 'RBR_data/20190927/';
files = {'20190927_LL_124127_data.txt','20190927_PN_124168_data.txt','20190927_SB_124163_data.txt','20190927_SL_124166_data.txt'};
labels = {'Lawsons Landing','Pelican Pt N','Seal Beach','Sacramento Landing'};

% Round Four
% prepend = 'RBR_data/20191124/';
% files = {'20191124_LL_124164_data.txt','20191124_SB_124165_data.txt','20191210_PN_124166_data.txt','20191210_TP_124162_data.txt'};
% labels = {'Lawsons Landing','Seal Beach','Pelican Point N','Tomasini Point'};

fs = 2; % 1/seconds
T = 1/fs;

rbr_times = {};
rbr_pressure = {};
rbr_depths_raw = {};
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
    sea_pressure = tmp{8};
    depth = tmp{9};
    
    timestamp = datetime(tmp{1},tmp{2},tmp{3},tmp{4},tmp{5},tmp{6});
    
    rbr_times{ii} = timestamp;
    rbr_pressures{ii} = pressure; % use end+1? 
    rbr_depths_raw{ii} = depth;
    
    fprintf('Done establishing times with RBR file %d of %d.\n',ii,length(files));
end

fprintf('Getting conditions...\n');
load conditions.mat

for kk = 1:length(files)
    fprintf('Adjusting depth data for sensor %d of %d.\n',kk,length(files));
    fixed_depths = zeros(size(rbr_pressures{kk}));
    for ll = 1:length(fixed_depths)
        progress_bar(ll,1,length(fixed_depths));
        [~,b_pressure_index] = min(abs(rbr_times{kk}(ll) - b_pressure.time));
        [~,temperature_index] = min(abs(rbr_times{kk}(ll) - temperature.time));
        [~,salinity_index] = min(abs(rbr_times{kk}(ll) - salinity.time));
        
        % rbr_pressures are direct pressure readings from the sensor. 
        
        g = 9.80665; % m/s^2
        
        % Note BML b_pressure is in mbar, RBR data comes in dbar
        % f_water_density takes pressure in bars
        % depth of water above sensor is (pressure - bpressure)/(density * g)
        
        rho = f_water_density(rbr_pressures{kk}(ll)/10, temperature.data(temperature_index), salinity.data(salinity_index));
        
        % We then convert pressures into Pascals for the final equation. 
        fixed_depths(ll) = (rbr_pressures{kk}(ll)*10000 - b_pressure.data(b_pressure_index)*100)/(g*rho);

    end
    
    rbr_depths_adjusted{kk} = fixed_depths;
    
end

save rbr_data_deployment_3.mat files labels rbr_times rbr_pressures rbr_depths_adjusted rbr_depths_raw

%% Extra Zone
% 
% prepend = 'RBR_data/20190927/';
% files = {'20190927_LL_124127.rsk','20190927_PN_124168.rsk','20190927_SB_124163.rsk','20190927_SL_124166.rsk'};
% labels = {'Lawsons Landing','Pelican Point North','Seal Beach','Sacramento Landing'};
% 
% rbr_times = {};
% rbr_depths = {};
% rbr_pressures = {};
% rbr_depths_adjusted = {};
% 
% for ii = 1:length(files)
%     rsk = RSKopen([prepend,files{ii}]);
%     rsk = RSKreaddata(rsk);
%     
%     rbr_times{ii} = datestr(rsk.data.tstamp);
%     rbr_depths{ii} = rsk.data.values;
%     progress_bar(ii,1,length(files));
% end
% 
% get_ocean_conditions

