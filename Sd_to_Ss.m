function Ss = Sd_to_Ss(ts,mab,W,Sd,varargin)
% SD_TO_SS adjusts a pressure timeseries (ts), collected by a pressure
% transducer some distance from the bed (mab) across a range of pulsations
% (W = 2*pi*f) to transform the power spectrum from the subsurface 
% timeseries (Sd) into a power spectrum that represents the surface (Ss). 
% The 4th argument argument is an optional value for N (typically 
% 0.90-1.07 or 1). All per Bishop & Donelan 1987, Coastal Engineering. 

% Lukas WinklerPrins
% lukas_wp@berkeley.edu
% 3 April 2020

% NEED TO FINISH HIGH-FREQUENCY CUTOFF STUFF

    if nargin == 5
        N = varargin{1};
    else
        N = 1;
    end
    % N is the empirical correction factor
    
    avg_mab = mean(ts);
    kh = qkhfs(W,avg_mab);
    kc = qkhfs(W,mab);
    Kp = cosh(kc)./cosh(kh);
    % Kp is the pressure response factor
    Ss = Sd.*(N./Kp).^2;
    
    % Find index of high-frequency cutoff
    g = 9.80665;
%     high_freq_cutoff = sqrt(g/(4*pi
end