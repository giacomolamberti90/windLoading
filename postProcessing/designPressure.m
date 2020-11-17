function [cp_design, area] = designPressure(str, tile, ang, window, POE, method)

    % Function that reads tile and POE 
    % returns cp design relative to POE, as a function of the area of the
    % panels

%%
if str == 'on'
    
    figure
        hold on
        scatter(tile.coords(:,1), tile.coords(:,2), 20, 'k', 'FaceColor', 'flat');

        xlabel('$$x[m]$$','interpreter','latex')
        ylabel('$$y[m]$$','interpreter','latex')
        pbaspect([1 1 1]); %axis([0 1 0 2])
        set(gca, 'Color', [0.8 0.8 0.8])
        set(gca, 'fontsize', 20)
end

delta = [2, 3, 4, 6, 9, 13, 18, 25];

for d = 1:length(delta)

    % panels sides
    dx = delta(d)*0.01;
    dy = delta(d)*0.015;

    % area of panel
    area(d) = dx*dy;

    % panel coordinates
    if ang == 180
        if tile.taps(1,1) == 'A'
            panel = [0, dx, 2.0-dy, 2.0];
        elseif tile.taps(1,1) == 'B'
            panel = [0, dx, 1.0-dy/2, 1.0+dy/2];
        end
    else
        if tile.taps(1,1) == 'A'
            panel = [1.0-dx, 1.0, 2.0-dy, 2.0];
        elseif tile.taps(1,1) == 'B'
            panel = [1.0-dx, 1.0, 1.0-dy/2, 1.0+dy/2];
        end
    end
    
    % compute area-averaged pressure
    tile.areaAverage = areaAveragedPressure(tile, panel);
    
    if method == 'cook'
        % compute design cp
        cp_design(d) = gumbel(tile, window, POE, [0.5,0.5,0.5]);
        %cp_design(d) = mean(tile.areaAverage) - 3 * std(tile.areaAverage);
        
    elseif method == 'safe'
        
        [mean_AA, var_AA] = areaAveragedStatistics(tile, panel);
        cp_design(d) = mean_AA - window * sqrt(var_AA);
    end
    % plot panel
    if str == 'on'
        rectangle('Position', [panel(1), panel(3), dx, dy], 'EdgeColor','r')
    end    
end
