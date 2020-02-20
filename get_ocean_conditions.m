% Lukas WinklerPrins, lukas_wp@berkeley.edu
% Last Modified 27 January 2020

%%%%% ATMOSPHERIC PRESSURE
fid = fopen('external_data/old/bml_barometric_pressure_2019_hourly.txt','rt');
P = textscan(fid,'%d,%d-%d-%d %d:%d:%d,%d,%f');
bpressure = P{9}; % in millibars, hourly averages
n = length(bpressure);
pressure_time = datetime(P{2},P{3},P{4},P{5},P{6},P{7});
fclose(fid);

%%%%% Seawater Density
fid = fopen('external_data/old/bml_seawater_density_2019_hourly.txt','rt');
P = textscan(fid,'%d,%d-%d-%d %d:%d:%d,%d,%f');
water_density = P{9} + 1000; % kg/m3
density_time = datetime(P{2},P{3},P{4},P{5},P{6},P{7});
fclose(fid);

%%%%% TIDES 
% Tide Data in GMT from Pt Reyes NOAA Buoy, in m, by minute
fid = fopen('external_data/old/9415020_tide.txt','rt');
P = textscan(fid,'%d/%d/%d %s %d:%d %f','HeaderLines',14);
tide = P{7}; % depths 

tide_time = datetime(P{1},P{2},P{3},P{5},P{6},zeros(size(P{6}))) - hours(7); % Seven Hour Adjustment, GMT to PDT
fclose(fid);

%%%% SWELL AND WIND
% Swell Data is only at 40min after the hour from 46013 (needs skipping)
% Swell Data is only every 30min after the hour from 46214 - won't use.
% Using NOAA Buoy 46013, hourly on the 40min (wvht ~= 99.00 to select actual measurements)
fid = fopen('external_data/old/46013_swell.txt','rt');
P = textscan(fid,'%d %d %d %d %d %d %f %f %f %f %f %d %f %f %f %f %f %f','HeaderLines',2);
wvht_temp = P{9}; %wvht = [wvht ~= 99.00];  
wdir_temp = P{12};  
dwpd_temp = P{10};
windsp_temp = P{7};
winddir_temp = P{6};
I = wvht_temp ~= 99.00;
wvht = [];
wdir = [];
dwpd = [];
windsp = [];
winddir = [];
temp_time = [];
for ii = 1:length(I)
    if I(ii)
        wvht = [wvht; wvht_temp(ii)];
        wdir = [wdir; wdir_temp(ii)];
        dwpd = [dwpd; dwpd_temp(ii)];
        windsp = [windsp; windsp_temp(ii)];
        winddir = [winddir; winddir_temp(ii)]; 
        temp_time = [temp_time; P{1}(ii),P{2}(ii),P{3}(ii),P{4}(ii),P{5}(ii),0];
    end
end

swell_time = datetime(temp_time);
fclose(fid);



