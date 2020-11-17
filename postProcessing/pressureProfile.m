function pressureProfile(tile, stat, label, col)

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
    axLim(2) = 0.5;
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
    axLim(2) = 0;    
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

% define rows of taps
if tile.taps(1,1) == 'A'
    for i = 13%1:15%13
        ind{i}  = find(str2num(tile.taps(:,2:3)) == i);
        
        %subplot(3,5,i)               
        if tile.CI == 'on' 
            if (stat == 'mean' | stat == 'stan' | stat == 'peak')
                hold on
                x  = [tile.coords(ind{i},1); flipud(tile.coords(ind{i},1))];
                y  = [tile.CI95(ind{i},1);   flipud(tile.CI95(ind{i},2))];                

                fill(x, y, col, 'linestyle', 'none') 
                alpha(0.5)
                %errorbar(tile.coords(ind{i},1), var(ind{i}),...
                %         var(ind{i})-tile.CI95(ind{i},1),...
                %         var(ind{i})+tile.CI95(ind{i},2), label, 'markersize', 10)
            end
        end
        hold on
        if label(1) == 'x'
            plot(tile.coords(ind{i},1), variable(ind{i}), label, 'markersize', 20)
        else
            plot(tile.coords(ind{i},1), variable(ind{i}), label)
        end
        
        %title(sprintf('y = %6.3f m', mean(tile.coords(ind{i},2))))    
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
    
elseif tile.taps(1,1) == 'B'
    
    for i = 8%1:19
        ind{i} = find(str2num(tile.taps(:,4:5)) == i);
  
        %subplot(4,5,i)
        if tile.CI == 'on'
            if (stat == 'mean' | stat == 'stan' | stat == 'peak')
                hold on
                x  = [tile.coords(ind{i},1); flipud(tile.coords(ind{i},1))];
                y  = [tile.CI95(ind{i},1);   flipud(tile.CI95(ind{i},2))];                

                fill(x, y, col, 'linestyle', 'none')           
                alpha(0.5)
                %errorbar(tile.coords(ind{i},1), var(ind{i}),...
                %         var(ind{i})-tile.CI95(ind{i},1),...
                %         var(ind{i})+tile.CI95(ind{i},2), label, 'markersize', 10)
            end
        end
        hold on
        if label(1) == 'x'
            plot(tile.coords(ind{i},1), variable(ind{i}), label, 'markersize', 20)
        else
            plot(tile.coords(ind{i},1), variable(ind{i}), label)
        end
        
        %title(sprintf('y = %6.3f m', mean(tile.coords(ind{i},2))))    
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
    
elseif tile.taps(1,1) == 'D'
    
    for i = 1:19%8
        ind{i} = find(str2num(tile.taps(:,4:5)) == i);
  
        subplot(4,5,i)
        if tile.CI == 'on'
            if (stat == 'mean' | stat == 'stan')            
                hold on
                x  = [tile.coords(ind{i},1); flipud(tile.coords(ind{i},1))];
                y  = [tile.CI95(ind{i},1);   flipud(tile.CI95(ind{i},2))];                

                fill(x, y, col, 'linestyle', 'none')           
                alpha(0.5)
                %errorbar(tile.coords(ind{i},1), var(ind{i}),...
                %         var(ind{i})-tile.CI95(ind{i},1),...
                %         var(ind{i})+tile.CI95(ind{i},2), label, 'markersize', 10)
            end
        end
        hold on
        plot(tile.coords(ind{i},1), variable(ind{i}), label)
        
        %title(sprintf('y = %6.3f m', mean(tile.coords(ind{i},2))))    
        xlabel('x [m]', 'interpreter','latex', 'fontsize', 22)
        ylabel('$$C_p$$','interpreter','latex', 'fontsize', 22)
        
        set(gca,'fontsize',18)
        axis([0.025 0.3 0 1.5])
        
        box off      
    end 
    
end
