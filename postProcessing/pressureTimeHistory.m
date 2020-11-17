function varargout = pressureTimeHistory(str, tile, col, threshold, xLim, taps, variable)

    % Function that data structure of tile A/B or D, taps names (optional)
    % and returns: time history plot

%% select taps
switch nargin
    case 3
        variable  = tile.areaAverage;
        xLim = [0, max(tile.time)];
        threshold = -100;
    case 4
        variable  = tile.areaAverage;
        xLim = [0, max(tile.time)];
    case 5
        t = [];
        variable = tile.areaAverage;
        leg = '';
    case 6
        [~, ~, index] = intersect(taps, tile.taps(:,1:5), 'rows', 'stable');

        t = find(tile.timeHistory(:,index) < threshold*ones(size(tile.timeHistory(:,index))));
        sprintf('%i peaks below %i', length(t), threshold)
        sprintf('Minimum C_p of %6.4f', min(tile.timeHistory(t,index)))
        
        variable = tile.timeHistory(:,index);                
        leg = taps;
    
end

if str=='on' % center on the peak
            
            % plot taps
            hold on
            plot(tile.time-xLim(1), variable, col, 'linewidth', 2)
            %plot(tile.time(t), var(t), '^r', 'markersize', 10, 'MarkerFaceColor', [1 0 0])            

            % properties
            xlabel('$$t [s]$$','Interpreter','latex'); ylabel('$$C_p$$','Interpreter','latex');
            xlim(xLim-xLim(1))
            %ylim([-5 1])
            set(gca, 'fontsize', 22)      

            % legend
            %legend(leg, 'Location', 'SouthWest')
            
elseif str=='no'
    
            % plot taps
            hold on
            plot(tile.time, variable, col)
            %plot(tile.time(t), variable(t), '^r', 'markersize', 10, 'MarkerFaceColor', [1 0 0])            

            % properties
            xlabel('$$t [s]$$','Interpreter','latex'); ylabel('$$C_p$$','Interpreter','latex');
            xlim(xLim)
            %ylim([-5 1])
            set(gca, 'fontsize', 22)      

            % legend
            %legend(leg, 'Location', 'NorthEast')
end

% optional outputs
varargout{1} = xLim;       
            