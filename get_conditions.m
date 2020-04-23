% Lukas WinklerPrins, lukas_wp@berkeley.edu
% Last Modified 23 April 2020

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

%% NOAA BUOY MET STUFF
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

%% NOAA BUOY SPECTRAL
% creates "WIND" and "SWELL" structs
% Swell Direction, Period, Height; Wind Direction, Wind Speed
buoy_file = '46013_2019_spectral.txt';
fid = fopen([path_prefix buoy_file],'rt');
P = textscan(fid, '%d %d %d %d %d %f %f %f %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',1);

buoy_spectra.file = buoy_file;
buoy_spectra.interval = 'hourly at 40min';
buoy_spectra.spectra = [];
temp_time = [];
buoy_spectra.freq = [.0200  .0325  .0375  .0425  .0475  .0525  .0575  .0625  .0675  .0725  .0775  .0825  .0875  .0925  .1000  .1100  .1200  .1300  .1400  .1500  .1600  .1700  .1800  .1900  .2000  .2100  .2200  .2300  .2400  .2500  .2600  .2700  .2800  .2900  .3000  .3100  .3200  .3300  .3400  .3500  .3650  .3850  .4050  .4250  .4450  .4650  .4850];

for ii = 1:length(P{1})
    progress_bar(ii,1,length(P{1}));
    temp = [];
    for jj = 6:52
        temp(jj-5) = P{jj}(ii);
    end
    buoy_spectra.spectra = [buoy_spectra.spectra; temp];
    temp_time = [temp_time; P{1}(ii),P{2}(ii),P{3}(ii),P{4}(ii),P{5}(ii),0];
end

buoy_spectra.time = datetime(temp_time);

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
