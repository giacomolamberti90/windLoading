function cp_peak = gumbel(tile, window, POE, npeaks, start)

    % Function that reads tile and number of windows and
    % returns gumbel distribution and POE% POE 
    
%% Minima -----------------------------------------------------------------
window_fs = 36; % 36 seconds in 1/50 scale and 1/3 velocity is 10 minutes in full scale   

W = floor(tile.time(end)/window); % 6 seconds in 1/50 scale is 2 minutes in full scale

% divide signal in N windows same number of samples
N = floor(length(tile.time)/W);

% reduced variate
% npeaks = W;

p = (1:npeaks)/(npeaks+1);
y = -log(-log(p));

k = 0;
for w = start + 1 : start + npeaks
    k = k + 1;
    % locate peaks
    peaks(k) = max(abs(tile.areaAverage(1+(w-1)*N : w*N)));
end
% fit Gumbel distribution
params = polyfit(y, sort(peaks), 1);

a = params(1); % slope
U = params(2); % intercept

U = U + a * log(window_fs/window); % correction Cook-Mayne 1980

% peak relative to POE 
cp_peak = - (U - a * log(-log(1-POE)));