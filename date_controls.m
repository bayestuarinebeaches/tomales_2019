%% Sensor Data & Dates

% Bring in labels{}, rbr_depths{}, rbr_times{}, rbr_pressures{}, rbr_depths_adjusted{}

%%% For Botany Bay Data
% start = [2018,6,22,12,0,0];
% start_time = datenum(start);
% end_time = datenum(2018,7,24,0,0,0);
% cmab = [0.3 0.15 0.15]; % RBR 1, 2, 5. RBR2's is guessed
% load bb_data_1-2-5.mat
% rbr_depths_adjusted = rbr_depths; % Because data not corrected. 
% 
% start_time = datenum(2018,7,13,0,0,0);
% end_time = datenum(2018,7,16,0,0,0);
    
%%% For FIRST Sensor Set
% NEED TO GET RBR DATA INTO MAT...
% start_time = datenum(2019,6,4,12,0,0);
% end_time = datenum(2019,7,16,0,0,0);
% %        LL  PN  PS  SB  WB  SL   TB
% cmab = [460, 20, 14, 47, 25, 116, 30];

%%% For SECOND Sensor Set
% load rbr_data_deployment_2.mat
% start = [2019,7,19,0,0,0];
% start_time = datenum(start);
% end_time = datenum(2019,8,29,0,0,0);
% % LL, PS, SL, TB, TP, WB
% % N.B. Wall Beach got mangled on Aug 5?try to go before or after this
% % Tomasini only deployed 8/1
% cmab = [460, 14, 116, 30, 3, 25];
% start_time = datenum(2019,8,6,0,0,0);

% Random period to examine
% start_time = datenum(2019,7,26,0,0,0);
% end_time = datenum(2019,7,28,0,0,0);

% Good sea breeze pattern days... 
% start_time = datenum(2019,8,17,0,0,0);
% end_time = datenum(2019,8,24,0,0,0);

% Visible Seiche (LL)...?
% start_time = datenum(2019,8,9,0,0,0);
% end_time = datenum(2019,8,10,0,0,0);

% Peak in Seiching noticeable again (LL?)
% start_time = datenum(2019,7,24,6,0,0);
% end_time = datenum(2019,7,25,18,0,0);

%%% For THIRD Sensor Set
load rbr_data_deployment_3.mat
start = [2019,8,30,0,0,0];
start_time = datenum(start);
end_time = datenum(2019,9,27,0,0,0);
% Note big time jump error in LL on Sept 23rd
% Note Tomasini Point data goes UNTIL 22 NOVEMBER? uses same data as 4th deploy
% LL, PPN, SB, SL, TP
cmab = [460, 20, 14, 116, 3];
end_time = datenum(2019,9,22,0,0,0);

% start_time = datenum(2019,9,4,22,30,0);
% end_time = datenum(2019,9,5,22,15,0);

% start_time = datenum(2019,9,27,12,0,0);
% end_time = datenum(2019,9,28,0,0,0);

% hard "start-up" of winds
% start_time = datenum(2019,9,24,12,0,0);
% end_time = datenum(2019,9,27,6,0,0);

%%% For FOURTH Sensor Set
% start = [2019,9,28,0,0,0];
% start_time = datenum(start);
% end_time = datenum(2019,11,21,0,0,0);
% % load rbr_data_deployment_4.mat
% load rbr_data_deployment_4_fixed.mat % Has fixed pressure records. 
% % LL, SB, PPN, TP
% % NOTE LAWSONS LANDING TIMESERIES HERE IS MESSED UP, CAME LOOSE SOME WAY IN...
% cmab = [460, 14, 20, 3];
% % end_time = datenum(2019,10,20,0,0,0);

%% OK Actually Now, from 4th Deployment

% Oct 27, Wind Day
% start_time = datenum(2019,10,27,8,0,0);
% end_time = datenum(2019,10,28,0,0,0);

% Nov 7, Swell Day
% start_time = datenum(2019,11,7,12,0,0);
% end_time = datenum(2019,11,8,12,0,0);

% Oct 18, "Combo" Day
% start_time = datenum(2019,10,18,0,0,0);
% end_time = datenum(2019,10,19,0,0,0);

% Consistently Windy over a couple tides
% start_time = datenum(2019,9,28,8,0,0);
% end_time = datenum(2019,9,29,8,0,0);
% 
% start_time = datenum(2019,9,27,0,0,0);
% end_time = datenum(2019,9,29,0,0,0);


% Other
% start_time = datenum(2019,9,24,8,0,0);
% end_time = datenum(2019,9,26,8,0,0);


%% Rubbish Bin

% Sept 28-29 - Big Wind Event, Spring Tide, No particular swell
% start_time = datenum(2019,9,28,12,0,0);
% end_time = datenum(2019,9,29,12,0,0);

% % Oct 17-Oct 19 - Windy, swell arrives, some big igw motions
% start_time = datenum(2019,10,17,0,0,0);
% end_time = datenum(2019,10,19,12,0,0);
% 
% start_time = datenum(2019,10,18,20,0,0);
% end_time = datenum(2019,10,19,0,0,0);

% start_time = datenum(2019,10,27,13,0,0);
% end_time = datenum(2019,10,27,17,0,0);

% % Nov 15-16 - Windy, bigger swell arrives, Lukas' Birthday
% start_time = datenum(2019,11,16,0,0,0);
% end_time = datenum(2019,11,16,18,0,0);

% % Nov 20-21 - bigger waves, from local wind? 
% start_time = datenum(2019,11,19,0,0,0);
% end_time = datenum(2019,11,20,12,0,0);

% Nov 7-8 - 12s consistent swell, very light winds
% start_time = datenum(2019,11,7,12,0,0);
% end_time = datenum(2019,11,8,12,0,0);

% N wind on Ebb tides
% start_time = datenum(2019,10,16,18,0,0);
% end_time = datenum(2019,10,18,6,0,0);

% S wind on Ebb tides
% start_time = datenum(2019,10,9,0,0,0);
% end_time = datenum(2019,10,11,18,0,0);

% N wind on Flood tides
% start_time = datenum(2019,10,7,0,0,0);
% end_time = datenum(2019,10,10,0,0,0);

% % S wind on Flood tides
% start_time = datenum(2019,10,1,0,0,0);
% end_time = datenum(2019,10,3,12,0,0);