% Lukas WinklerPrins
% lukas_wp@berkeley.edu

% Last Edited 14 January 2020

prepend = 'RBR_data/20190927/';
files = {'20190927_LL_124127.rsk','20190927_PN_124168.rsk','20190927_SB_124163.rsk','20190927_SL_124166.rsk'};
labels = {'Lawsons Landing','Pelican Point North','Seal Beach','Sacramento Landing'};

rbr_times = {};
rbr_depths = {};
rbr_pressures = {};
rbr_depths_adjusted = {};

for ii = 1:length(files)
    rsk = RSKopen([prepend,files{ii}]);
    rsk = RSKreaddata(rsk);
    
    rbr_times{ii} = datestr(rsk.data.tstamp);
    rbr_depths{ii} = rsk.data.values;
    progress_bar(ii,1,length(files));
end

get_ocean_conditions

