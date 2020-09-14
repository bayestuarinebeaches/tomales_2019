%% Sensor Choice

% Remove this if using META_RUN.M
sensor_choice = 5;

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
make_spectra_plot = 1;
visualize_conditions = 1;
use_residual_spectra = 0;
do_energy_over_time = 1;
visualize_big_spectra = 1;
include_seiche = 0;

% Spectrum Choice: 0 for none
%                  1 for High/Trans/Low Tide
%                  2 for Morning/Midday/Evening
%                  3 for Flooding/Ebbing/Slack
%                  4 for High / Low / Ebbing / Flooding
spectra_choice = 4;

%% Timing Controls
% ALL THREE OF THE NEXT VARIABLES ARE IN HOURS
ea_spacing = 1; %6; % What is the spacing between each ensemble?
window_length = 3; %3 %24; % Each ensemble represents how many hours? 
instance_length = 3/4; %3; % What is the length of each instance in an ensemble?

% 1 /  3 / 0.75 has settled out to be the default
% 6 / 18 / 6 is useful for resolving seiche hump w/o tidal consideration

% 24 / 24 / 0.75 ???

