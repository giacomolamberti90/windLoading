function varargout = pressurePDF(tile, col, taps)

    % Function that data structure of tile A/B or D, taps names (optional)
    % and returns: time history plot

%% select taps
switch nargin
    case 2
        var    = tile.areaAverage;
        xi     = linspace(min(var), max(var), 100);
        [f, x] = ksdensity(var, xi);     
    case 3
        [~, ~, index] = intersect(taps, tile.taps(:,1:5), 'rows', 'stable');
        
        var    = tile.timeHistory(:,index);
        xi     = linspace(-2, 0, 100);
        [f, x] = ksdensity(var, xi);            
end
% skewness
disp(sprintf('skewness = %f', skewness(var)))

% plot pdf
plot(x, f, 'linewidth', 1.5, 'color', col)
%histogram(var, 'Normalization', 'pdf', 'FaceColor', col, 'EdgeColor', 'none')

% properties
ylabel('$$f$$','Interpreter','latex'); xlabel('$$C_{p,AA}$$','Interpreter','latex');
set(gca, 'fontsize', 18)       