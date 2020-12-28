%% Sensor Choice

% Remove this if using META_RUN.M
sensor_choice = input('Which sensor do you want? ');

%% Mechanical Options
% "1" turns them on, "0" turns them off
adjust_pressure = 1; % Default is 1
variance_preserving = 1; % Default is 1
use_t_tide = 0; % Default is 1, turn off if very short time window...
use_windowing = 1; % Default is 0
use_median_power = 0; % Not functional yet (keep at 0).
low_pass_filter = 0; % Default is 0
pre_smooth = 0; % Default is 0

%% Visualization Options
% "1" turns them on, "0" turns them off
make_contour_graph = 1;
make_conditions_plot = 1;
make_spectral_plot = 1;
make_wave_height_plot = 1;
make_energy_bar_chart = 0;
make_polar_plot = 0;
make_chop_plot = 0;
make_swell_plot = 0;
make_wind_start_plot = 0;
make_energy_flux_contour_graph = 0;
plot_only_mean_spectrum = 0;

use_residual_spectra = 0; % Default is 0
include_seiche = 0; % Default is 0

%% Timing Controls
% ALL THREE OF THE NEXT VARIABLES ARE IN HOURS
ea_spacing = 3; %6; % What is the spacing between each ensemble?
window_length = 6; %3 %24; % Each ensemble represents how many hours? 
instance_length = 0.75; %3; % What is the length of each instance in an ensemble?

% 1 /  3 / 0.75 has settled out to be the default!!! 
% 6 / 18 / 6 is useful for resolving seiche hump w/o tidal consideration

% 24 / 24 / 0.75 ???

%% Extra Controls

window_slope_cutoff = 0.15; % 0.15 seems to strike a good balance between the 4 categories. In m/hour.

wind_option = 5;
% 1 uses BOON, 2 uses TBB, 3 uses HIOC, 4 uses TBOC
% 5 uses TBB for speed and BOON for wind direction! 