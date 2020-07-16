% Lukas WinklerPrins
% lukas_wp@berkeley.edu
% Last edited 15 July 2020

% Lawsons Landing
% 2 Jun - 17 Jul
% 17 Jul - 29 Aug
% 29 Aug - 27 Sep with big time error event on Sept 23rd
% 27 Sep - 24 Nov sensor came loose somewhere in the deployment (~Oct 20th?)
% 24 Nov - 6 Feb
LL = [datetime(2019,6,2) datetime(2019,7,17);
    datetime(2019,7,17) datetime(2019,8,29);
    datetime(2019,8,29) datetime(2019,9,23);
    datetime(2019,9,23) datetime(2019,9,27);
    datetime(2019,9,27) datetime(2019,10,20);
    datetime(2019,10,20) datetime(2019,11,24);
    datetime(2019,11,24) datetime(2020,2,6)];

LL_style = {'ks-'; 'ks-'; 'ks-'; 'k:'; 'ks-'; 'k:'; 'ks-'};

% Seal Beach
% 2 Jun - 17 Jul
% x
% 29 Aug - 27 Sep
% 27 Sep - 24 Nov
% 24 Nov - 6 Feb
% 6 Feb - 11 May
SB = [datetime(2019,6,2) datetime(2019,7,17);
    datetime(2019,8,29) datetime(2019,9,27);
    datetime(2019,9,27) datetime(2019,11,24);
    datetime(2019,11,24) datetime(2020,2,6)];

SB_style = {'ks-'; 'ks-'; 'ks-'; 'ks-'};

% Wall Beach
% 2 Jun - 17 Jul
% 17 Jul - 29 Aug with "mangling" event on Aug 5th (data are OK before & after, but spike day-of)
WB = [datetime(2019,6,2) datetime(2019,7,17);
    datetime(2019,7,17) datetime(2019,8,4);
    datetime(2019,8,4) datetime(2019,8,6);
    datetime(2019,8,6) datetime(2019,8,29)];
WB_style = {'ks-'; 'ks-'; 'k:'; 'ks-'};

% Pelican North
% 2 Jun - 17 Jul
% 17 Jul - 29 Aug
% 29 Aug - 27 Sep
% 27 Sep - 24 Nov
% 24 Nov - 6 Feb
% 6 Feb - 11 May
PN = [datetime(2019,6,2) datetime(2019,7,17);
    datetime(2019,7,17) datetime(2019,8,29);
    datetime(2019,8,29) datetime(2019,9,27);
    datetime(2019,9,27) datetime(2019,11,24);
    datetime(2019,11,24) datetime(2020,2,6)];
PN_style = {'ks-'; 'ks-'; 'ks-'; 'ks-'; 'ks-'};

% Pelican South
% 2 Jun - 17 Jul
% 17 Jul - 29 Aug
PS = [datetime(2019,6,2) datetime(2019,7,17);
    datetime(2019,7,17) datetime(2019,8,29)];
PS_style = {'ks-'; 'ks-'};

% Tomales Beach
% 2 Jun - 17 Jul
% 17 Jul - 29 Aug
TB = [datetime(2019,6,2) datetime(2019,7,17);
    datetime(2019,7,17) datetime(2019,8,19)];
TB_style = {'ks-'; 'ks-'};

% Sac Landing
% 2 Jun - 17 Jul
% 17 Jul - 29 Aug
% 29 Aug - 27 Sep
% 27 Sep - 24 Nov
% 24 Nov - 6 Feb
SL = [datetime(2019,6,2) datetime(2019,7,17);
    datetime(2019,7,17) datetime(2019,8,19);
    datetime(2019,8,29) datetime(2019,9,27);
    datetime(2019,9,27) datetime(2019,11,24);
    datetime(2019,11,24) datetime(2020,2,6)];
SL_style = {'ks-'; 'ks-'; 'ks-'; 'ks-'; 'ks-'};

% Tomasini Point
% 1 Aug - 29 Aug
% 29 Aug - 27 Sep
% x
% 24 Nov - 6 Feb
% 6 Feb - 11 May
TP = [datetime(2019,8,1) datetime(2019,8,29);
    datetime(2019,8,29) datetime(2019,9,27);
    datetime(2019,11,24) datetime(2020,2,6)];
TP_style = {'ks-'; 'ks-'; 'ks-'};

sensors = {LL,SB,WB,PN,PS,TB,SL,TP};
styles = {LL_style,SB_style,WB_style,PN_style,PS_style,TB_style,SL_style,TP_style};
labels = {'Lawsons Landing','Seal Beach','Wall Beach','Pelican Point N','Pelican Point S','Tomales Beach','Sacramento Landing','Tomasini Point'};

figure
hold on
for jj = 1:length(sensors)
    for ii = 1:size(sensors{jj},1)
        plot(sensors{jj}(ii,:),[(length(sensors)-jj+1) (length(sensors)-jj+1)],styles{jj}{ii});
    end
    text(datetime(2019,6,15),(length(sensors)-jj+1.3),labels{jj});
end

ylim([0 9]);