path_prefix = 'external_data/';
station_file = '2019_HIOCweather.txt';
fid = fopen([path_prefix station_file],'rt');
M = textscan(fid, '%d/%d/%d %d:%d %s %s %s %s %s %f %s %f %f %s %s %s %s %f %f %f %s %s %f %f %f %f %f %f %f %d %f %f','HeaderLines',2);

hioc_wind.file = station_file;
hioc_wind.interval = 'Every 30min on the hour';
hioc_wind.spd = [];
hioc_wind.spd_unit = 'mph'; % I think? 
hioc_wind.dir = [];
hioc_wind.dir_unit = 'deg';
temp_time = [];

temp = M{12};
I = ~strcmp(M{12},'---'); % Only using points with direction data
for ii = 1:length(I)
    if I(ii)
        hioc_wind.spd = [hioc_wind.spd; M{11}(ii)];
        switch char(M{12}(ii))
            case 'N'
                hioc_wind.dir = [hioc_wind.dir; 360];
            case 'NNE'
                hioc_wind.dir = [hioc_wind.dir; 22.5];
            case 'NE'
                hioc_wind.dir = [hioc_wind.dir; 45];
            case 'ENE'
                hioc_wind.dir = [hioc_wind.dir; 67.5];
            case 'E'
                hioc_wind.dir = [hioc_wind.dir; 90];
            case 'ESE'
                hioc_wind.dir = [hioc_wind.dir; 112.5];
            case 'SE'
                hioc_wind.dir = [hioc_wind.dir; 135];
            case 'SSE'
                hioc_wind.dir = [hioc_wind.dir; 157.5];
            case 'S'
                hioc_wind.dir = [hioc_wind.dir; 180];
            case 'SSW'
                hioc_wind.dir = [hioc_wind.dir; 202.5];
            case 'SW'
                hioc_wind.dir = [hioc_wind.dir; 225];
            case 'WSW'
                hioc_wind.dir = [hioc_wind.dir; 247.5];
            case 'W'
                hioc_wind.dir = [hioc_wind.dir; 270];
            case 'WNW'
                hioc_wind.dir = [hioc_wind.dir; 292.5];
            case 'NW'
                hioc_wind.dir = [hioc_wind.dir; 315];
            case 'NNW'
                hioc_wind.dir = [hioc_wind.dir; 337.5];
        end
        temp_time = [temp_time; (M{3}(ii)+2000) M{1}(ii) M{2}(ii) M{4}(ii) M{5}(ii) 0];
    end
end

hioc_wind.time = datetime(temp_time);

save hioc_wind.mat hioc_wind
        