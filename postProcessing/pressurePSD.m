function varargout = pressurePSD(tile, fsamp, str, col, taps, variable)

    % Function that data structure of tile A/B or D, sampling frequency,
    % taps names and returns: pressure power spectrum

%% select taps
switch nargin
    case 4
        variable   = (tile.areaAverage - mean(tile.areaAverage));
        sigma = std(tile.areaAverage);
    case 5
        [~, ~, index] = intersect(taps, tile.taps(:,1:5), 'rows', 'stable');        
        variable   = (tile.timeHistory(:,index)-tile.mean(index));
        sigma = tile.std(index);
end

%% raw spectrum
% Nsamp = length(tile.time);
% f     = linspace(1/tile.time(end), fsamp/2, 1+Nsamp/2);
% df    = f(2)-f(1);
% pk_f  = fft(tile.timeHistory(:,index))/Nsamp;
%  
% if (mod(Nsamp,2) == 0) % N even
%     for i = 1:Nsamp/2
%        S(i) = 2*abs(pk_f(i)).^2/df;
%     end
%     % Nyquist frequency
%     S(Nsamp/2+1) = abs(pk_f(Nsamp/2+1)).^2/df;
% else
%     for i = 1:(Nsamp-1)/2+1
%        S(i) = 2*abs(pk_f(i)).^2/df;
%     end  
% end
%  
% % plot taps
% loglog(f(2:end), filter([0.25 0.25 0.25 0.25], 1, S(2:end)), 'linewidth', 1.5)
% 
% % properties
% ylabel('$$S[]$$','Interpreter','latex'); xlabel('$$f[Hz]$$','Interpreter','latex');
% set(gca, 'fontsize', 18)      
% 
% return

%% pwelch
[S, f] = pwelch(variable, [], [], [], fsamp);

% plot taps
if str=='dim'
    loglog(f, S, 'linewidth', 1.5, 'color', col)

    % properties
    xlabel('$$f [Hz]$$','Interpreter','latex'); 
    ylabel('$$S [Pa^2s]$$','Interpreter','latex');
    set(gca, 'fontsize', 18)
    box off 
    
elseif str=='adi'
    
    loglog(f*0.3/tile.U, S.*f/sigma^2, 'linewidth', 1.5, 'color', col)    
    ylim([1e-7 10])
    
    % properties
    ylabel('$$fS/C_{p,rms}$$','Interpreter','latex'); 
    xlabel('$$fB/U$$','Interpreter','latex');
    set(gca, 'fontsize', 18)
    box off

elseif str=='vor'  
    
    [S, f] = pwelch(variable, 1000, [], [], fsamp);
    
    S_ad = S.*f/sigma^2;
    f_ad = f*0.3/tile.U;
    
    strouhal = f_ad(find(S_ad == max(S_ad)));
    
    disp(sprintf('Strouhal = %6.3f', strouhal));
    disp(sprintf('frequency of vortex shedding = %6.3f', strouhal*tile.U/0.3));
    
    hold on
    plot(f_ad, S_ad, 'linewidth', 1.5, 'color', col)
    plot(strouhal, S_ad(find(S_ad == max(S_ad))), '.r', 'markersize', 20)

    % properties
    ylabel('$$fS/\overline{p''^2}$$','Interpreter','latex'); 
    xlabel('$$fB/U$$','Interpreter','latex');
    set(gca, 'fontsize', 18)
    box off    
    
    xlim([0 20*strouhal])
end
