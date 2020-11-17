function pressureProfile_all(tile, str, stat, label, col)

    % Function that reads data structure of tile A/B or D 4 characters string 
    % containing the statistic we want to plot, the axis limits (optional)
    % and returns:
        
        %      axLim: axis limits (optional)

% select stastitic
if  stat == 'mean'
    if tile.CI == 'on'
        tile.CI95 = tile.CI95_mean;
    end
    variable = tile.mean;     
    axLim(1) = -2;
    axLim(2) = 0;
elseif stat == 'stan'
    if tile.CI == 'on'
        tile.CI95 = tile.CI95_std;
    end
    variable = tile.std;    
    axLim(1) = 0;
    axLim(2) = 0.3;
elseif stat == 'peak'
    if tile.CI == 'on'
        tile.CI95 = tile.CI95_peak;
    end       
    variable = tile.peak; 
    axLim(1) = -3;
    axLim(2) = -0.5;    
elseif stat == 'sob1'
    variable = tile.sobol(:,1); 
    axLim(1) = 0;
    axLim(2) = 1;
elseif stat == 'sob2'
    variable = tile.sobol(:,2); 
    axLim(1) = 0;
    axLim(2) = 1;
elseif stat == 'sob3'
    variable = tile.sobol(:,3); 
    axLim(1) = 0;
    axLim(2) = 1;
end

if str == 'A'
            
    ind = find(round(tile.coords(:,2),3) == 1.761); % A13XX
    %ind = find(round(tile.coords(:,2),3) == 1.947);
    [x_coord, xi] = sort(tile.coords(ind,1));

    if tile.CI == 'on' 
        if (stat == 'mean' | stat == 'stan' | stat == 'peak')
            hold on            
            CI_up = tile.CI95(ind,2);
            CI_low = tile.CI95(ind,1);
            
            x  = [x_coord; flip(x_coord)];
            y  = [CI_up(xi); flip(CI_low(xi))];

            fill(x, y, col, 'linestyle', 'none')           
            alpha(0.5)
        end
    end
    hold on
    plot(x_coord, variable(ind(xi)), label)

    xlabel('x [m]', 'interpreter','latex', 'fontsize', 22)
    if stat == 'mean'
        ylabel('$$C_p$$','interpreter','latex', 'fontsize', 22)
    elseif stat == 'stan'
        ylabel('$$C''_p$$','interpreter','latex', 'fontsize', 22)
    elseif stat == 'peak'
        ylabel('$$\check{C}_p$$','interpreter','latex', 'fontsize', 22)       
    end

    set(gca,'fontsize',18)
    axis([0.025 1 axLim(1) axLim(2)])

    box off        
    
elseif str == 'B'
            
    ind = find(round(tile.coords(:,2),3) == 0.998); % BXX10

    if tile.CI == 'on' 
        if (stat == 'mean' | stat == 'stan' | stat == 'peak')
            hold on
            [x, xi] = sort(tile.coords(ind,1));
            
            CI_up = tile.CI95(ind,2);
            CI_low = tile.CI95(ind,1);
            
            x  = [x; flip(x)];
            y  = [CI_up(xi); flip(CI_low(xi))];             

            fill(x, y, col, 'linestyle', 'none')           
            alpha(0.5)
        end
    end
    hold on
    plot(tile.coords(ind,1), variable(ind), label)

    xlabel('x [m]', 'interpreter','latex', 'fontsize', 22)
    if stat == 'mean'
        ylabel('$$C_p$$','interpreter','latex', 'fontsize', 22)
    elseif stat == 'stan'
        ylabel('$$C''_p$$','interpreter','latex', 'fontsize', 22)
    elseif stat == 'peak'
        ylabel('$$\check{C}_p$$','interpreter','latex', 'fontsize', 22)       
    end

    set(gca,'fontsize',18)
    axis([0.025 1 axLim(1) axLim(2)])

    box off  
    
 end