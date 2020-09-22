path_prefix = 'external_data/';
station_file = '2019_TBOCweather.txt';
fid = fopen([path_prefix station_file],'rt');
M = textscan(fid, '%d/%d/%d %d:%d %s %s %s %s %s %f %s %f %f %s %s %s %s %f %f %f %s %s %f %f %f %f %f %f %f %d %f %f','HeaderLines',2);

tboc_wind.file = station_file;
tboc_wind.interval = 'Every 30min on the hour';
tboc_wind.spd = [];
tboc_wind.spd_unit = 'mph'; % I think? 
tboc_wind.dir = [];
tboc_wind.dir_unit = 'deg';
temp_time = [];

temp = M{12};
I = ~strcmp(M{12},'---'); % Only using points with direction data
for ii = 1:length(I)
    if I(ii)
        tboc_wind.spd = [tboc_wind.spd; M{11}(ii)];
        switch char(M{12}(ii))
            case 'N'
                tboc_wind.dir = [tboc_wind.dir; 0];
            case 'NNE'
                tboc_wind.dir = [tboc_wind.dir; 22.5];
            case 'NE'
                tboc_wind.dir = [tboc_wind.dir; 45];
            case 'ENE'
                tboc_wind.dir = [tboc_wind.dir; 67.5];
            case 'E'
                tboc_wind.dir = [tboc_wind.dir; 90];
            case 'ESE'
                tboc_wind.dir = [tboc_wind.dir; 112.5];
            case 'SE'
                tboc_wind.dir = [tboc_wind.dir; 135];
            case 'SSE'
                tboc_wind.dir = [tboc_wind.dir; 157.5];
            case 'S'
                tboc_wind.dir = [tboc_wind.dir; 180];
            case 'SSW'
                tboc_wind.dir = [tboc_wind.dir; 202.5];
            case 'SW'
                tboc_wind.dir = [tboc_wind.dir; 225];
            case 'WSW'
                tboc_wind.dir = [tboc_wind.dir; 247.5];
            case 'W'
                tboc_wind.dir = [tboc_wind.dir; 270];
            case 'WNW'
                tboc_wind.dir = [tboc_wind.dir; 292.5];
            case 'NW'
                tboc_wind.dir = [tboc_wind.dir; 315];
            case 'NNW'
                tboc_wind.dir = [tboc_wind.dir; 337.5];
        end
        temp_time = [temp_time; (M{3}(ii)+2000) M{1}(ii) M{2}(ii) M{4}(ii) M{5}(ii) 0];
    end
end

tboc_wind.time = datetime(temp_time);

save tboc_wind.mat tboc_wind
        