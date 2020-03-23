function Ss = Sp_to_Ss(ts,mab,W,Sp,varargin)
% SP_TO_SS adjusts a pressure timeseries (ts), collected by a pressure
% transducer some distance from the bed (mab) across a range of pulsations
% (W = 2*pi*f) to transform the power spectrum from the subsurface 
% timeseries (Sp) into a power spectrum that represents the surface (Ss). 
% The 4th argument argument is an optional value for N (typically 
% 0.90-1.07 or 1). All per Bishop & Donelan 1987, Coastal Engineering. 

% Lukas WinklerPrins
% lukas_wp@berkeley.edu
% 23 February 2020

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
    Ss = Sp.*(N./Kp).^2;
end