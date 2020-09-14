function Ss = Sd_to_Ss(ts,mab,W,Sd,varargin)
% SD_TO_SS adjusts a pressure timeseries (ts), collected by a pressure
% transducer some distance from the bed (mab) across a range of pulsations
% (W = 2*pi*f) to transform the power spectrum from the subsurface 
% timeseries (Sd) into a power spectrum that represents the surface (Ss). 
% The 4th argument argument is an optional value for N (typically 
% 0.90-1.07 or 1). All per Bishop & Donelan 1987, Coastal Engineering. 

% Lukas WinklerPrins
% lukas_wp@berkeley.edu
% 29 April 2020

    if nargin == 5
        N = varargin{1};
    else
        N = 1;
    end
    % N is the empirical correction factor
    
    avg_depth = mean(ts);
    kh = qkhfs(W,avg_depth); % k * water depth
    kc = qkhfs(W,mab); % k * height above bed
    K_p = cosh(kc)./cosh(kh); % Pressure Response Factor
    Ss = Sd.*(N./K_p).^2;

end