% Lukas WinklerPrins, lukas_wp@berkeley.edu
% Last Modified 6 March 2020

path_prefix = 'external_data/';

% This is mostly using information that lives in a "metadata" spreadsheet.

%% BAROMETRIC PRESSURE
% creates "B_PRESSURE" struct
b_pressure.file = 'bml_barometric_pressure_2019_hourly.csv';
series_length = '8817';
b_pressure.time = readmatrix([path_prefix b_pressure.file],'Range',['B1:B' series_length],'OutputType','datetime');
b_pressure.data = readmatrix([path_prefix b_pressure.file],'Range',['D1:D' series_length],'OutputType','double');
b_pressure.unit = 'mbar';
b_pressure.interval = 'hourly';

%% SALINITY
% creates "SALINITY" struct
salinity.file = 'tbb_seawater_salinity_2019_hourly.csv';
series_length = '8343';
salinity.time = readmatrix([path_prefix salinity.file],'Range',['B1:B' series_length],'OutputType','datetime');
salinity.data = readmatrix([path_prefix salinity.file],'Range',['D1:D' series_length],'OutputType','double');
salinity.unit = 'psu';
salinity.interval = 'hourly';

%% TURBIDITY
turbidity.file = 'tbb_seawater_turbidity_2019_hourly.csv';
series_length = '8434';
turbidity.time = readmatrix([path_prefix turbidity.file],'Range',['B1:B' series_length],'OutputType','datetime');
turbidity.data = readmatrix([path_prefix turbidity.file],'Range',['D1:D' series_length],'OutputType','double');
turbidity.unit = 'NTU';
turbidity.interval = 'hourly';

%% WATER TEMPERATURE

temperature.file = 'tbb_seawater_temperature_2019_hourly.csv';
series_length = '8343';
temperature.time = readmatrix([path_prefix temperature.file],'Range',['B1:B' series_length],'OutputType','datetime');
temperature.data = readmatrix([path_prefix temperature.file],'Range',['D1:D' series_length],'OutputType','double');
temperature.unit = 'deg C';
temperature.interval = 'hourly';

%% NOAA BUOY STUFF
% creates "WIND" and "SWELL" structs
% Swell Direction, Period, Height; Wind Direction, Wind Speed
buoy_file = '46013_2019_stdmet.txt';
fid = fopen([path_prefix buoy_file],'rt');
P = textscan(fid, '%d %d %d %d %d %d %f %f %f %f %f %d %f %f %f %f %f %f','HeaderLines',2);

swell.file = buoy_file;
wind.file = buoy_file;
swell.interval = 'hourly at 40min';
wind.interval = 'hourly at 40min';
swell.hgt = []; % Significant Wave Height
swell.hgt_unit = 'meters';
swell.dir = []; % Direction
swell.dir_unit = 'deg';
swell.per = []; % Period
swell.per_unit = 'seconds';
wind.spd = []; % Wind Speed
wind.spd_unit = 'm/s';
wind.dir = []; % Direction
wind.dir_unit = 'deg';
temp_time = [];

temp = P{9};
I = temp ~= 99.00; % Only selecting the points with data
for ii = 1:length(I)
    if I(ii)
        swell.hgt = [swell.hgt; P{9}(ii)];
        swell.dir = [swell.dir; P{12}(ii)];
        swell.per = [swell.per; P{10}(ii)];
        wind.spd = [wind.spd; P{7}(ii)];
        wind.dir = [wind.dir; P{6}(ii)];
        temp_time = [temp_time; P{1}(ii),P{2}(ii),P{3}(ii),P{4}(ii),P{5}(ii),0];
    end
end

swell.time = datetime(temp_time);
wind.time = datetime(temp_time);

%% TIDES
% Not really crucial, but I guess we can have it... 

tide_lengths = [7441 6721 7441 7201 7441 7201 7441 7441 7201 7441 7201 7441];
    % the above INCLUDES counting the header line
    
tides.data = [];
tides.time = [];

for ii = 1:12
    tides.files{ii} = ['CO-OPS_9415020_wl_2019_' num2str(ii) '.csv'];
    temp_date = readmatrix([path_prefix tides.files{ii}],'Range',['A2:A' num2str(tide_lengths(ii))],'OutputType','char');
    temp_time = readmatrix([path_prefix tides.files{ii}],'Range',['B2:B' num2str(tide_lengths(ii))],'OutputType','char');
    temp_datetime = [];
    for jj = 1:length(temp_time)
        temp_datetime = [temp_datetime; datetime([temp_date{jj} ' ' temp_time{jj}],'InputFormat','yyyy/MM/dd HH:mm')];
    end
    tides.time = [tides.time; temp_datetime];
    tides.data = [tides.data; readmatrix([path_prefix tides.files{ii}],'Range',['E2:E' num2str(tide_lengths(ii))],'OutputType','double')];
end
tides.unit = 'meters';
tides.interval = '6 minutes';

save conditions.mat b_pressure salinity turbidity temperature swell wind tides 

%% TO-DO

% - Automate timeseries length acquisition
% - Automate grabbing data direct from web
