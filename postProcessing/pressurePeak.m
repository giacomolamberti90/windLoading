function cp_peak = pressurePeak(tile, window, POE, method)

    % Function that reads tile and number of windows and
    % returns gumbel distribution and POE% POE 

% Cook-Mayne
if method == 'cook'
    
    W_fs = 6; % 36 seconds in 1/50 scale amd 1/3 velocity is 10 minutes in full scale   

    W = floor(tile.time(end)/6); % 6 seconds in 1/50 scale is 2 minutes in full scale

    % divide signal in N windows same number of samples
    N = floor(length(tile.time)/W);

    % reduced variate
    npeaks = 10;
    
    p = (1:npeaks)/(npeaks+1);
    y = -log(-log(p));
    
    for tap = 1:size(tile.timeHistory, 2)
        k = 0;
        for w = 1 + (window-1) * npeaks : window * npeaks
            k = k + 1;
            % locate peaks
            peaks(k) = max(abs(tile.timeHistory(1+(w-1)*N:w*N, tap)));
        end
        % fit Gumbel distribution
        params = polyfit(y, sort(peaks), 1);
        
        a = params(1); % slope
        U = params(2); % intercept
        
        %U = U + a * log(W_fs/6); % correction Cook-Mayne 1980

        % peak relative to POE 
        cp_peak(tap) = - (U - a * log(-log(1-POE)));
    end

% safety factor
elseif method == 'safe'
    
    cp_peak = tile.mean' - window * tile.std';
    
elseif method == 'dave'
    
    for tap = 1:size(tile.timeHistory, 2)
        
        % spectral moments
        [S, f] = pwelch(tile.timeHistory(:,tap), [], [], [], 1/tile.time(2));
        
        m0 = trapz(f, S);
        m2 = trapz(f, f.^2.*S);
        
        % m0 = tfsmoment(tile.timeHistory(:,tap), tile.time, 0);
        % m2 = tfsmoment(tile.timeHistory(:,tap), tile.time, 2);
        
        nu = sqrt(m2/m0);
        T  = tile.time(end);
        
        g(tap) = sqrt(2 * log(nu*T)) + exp(1)/sqrt(2 * log(nu*T));
        
        cp_peak(tap) = tile.mean(tap) - g(tap) * tile.std(tap);
    end
    mean(g)

elseif method == 'wors'
    
    cp_peak = min(tile.timeHistory);
    
elseif method == 'ense'
    
    W = 5;
    N = floor(length(tile.time)/W);
    for tap = 1:size(tile.coords,1)
        for w = 1:W           
            % locate peaks
            peaks(w) = min(tile.timeHistory(1+(w-1)*N:w*N, tap));
        end

        % peak relative to POE 
        cp_peak(tap) = mean(peaks);
    end    
    
end
