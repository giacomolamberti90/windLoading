function varargout = pressureContour(str, tile, stat, axLim)

    % Function that reads string 'on'/'off' to plot or not taps position
    % and names, data structure of tile A/B or D 4 characters string 
    % containing the statistic we want to plot, the axis limits (optional)
    % and returns:
        
        %      axLim: axis limits (optional)
        
%% plot taps position        
if str=='on'
    
        figure('visible', 'on');
            hold on
            scatter(tile.coords(:,1), tile.coords(:,2), 20, 'k', 'FaceColor', 'flat');
            
            text(double(tile.coords(:,1)), double(tile.coords(:,2)), tile.taps, ...
                 'VerticalAlignment', 'middle', 'Rotation', 45, 'Margin', 15)
             
            xlabel('$$x[m]$$','interpreter','latex')
            ylabel('$$y[m]$$','interpreter','latex')
            set(gca, 'fontsize', 20)   
            
            fig2=gca;
            
            varargout{3} = fig2;
end

% select statistic        
if  stat == 'mean'
    var = tile.mean; 
elseif stat == 'stan'
    var = tile.std;
elseif stat == 'peak'
    var = tile.peak';    
elseif stat == 'sob1'
    var = tile.sobol(:,1);   
elseif stat == 'sob2'
    var = tile.sobol(:,2);    
elseif stat == 'sob3'
    var = tile.sobol(:,3);
elseif stat == 'mdel'
    var = tile.CI95_mean(:,2) - tile.CI95_mean(:,1);
elseif stat == 'sdel'
    var = tile.CI95_std(:,2) - tile.CI95_std(:,1);
elseif stat == 'spea'
    var = tile.CI95_peak(:,2) - tile.CI95_peak(:,1);
elseif stat == 'spod'
    var = real(tile.spod(:,1));
end

switch nargin
    case 3
        % axis limit
        axLim(1) = min(var);
        axLim(2) = max(var); 
end

% Triangulation
DT = delaunayTriangulation(double(tile.coords));

x = tile.coords(:,1);
y = tile.coords(:,2);

[X, Y] = meshgrid(linspace(min(x),max(x)), linspace(min(y),max(y)));

Z = griddata(double(x), double(y), double(var), double(X), double(Y));

        hold on
            % contour plot
            patch('Faces', DT.ConnectivityList, 'Vertices', [tile.coords, var],...
                  'FaceVertexCData', var, 'EdgeColor', 'none');
            
            %[C, h] = contourf(X, Y, Z, 7, ':');
            %h.LevelList = round(h.LevelList, 2);
            %clabel(C, h, 'FontWeight', 'bold');
            
            % plot taps
            plot3(tile.coords(:,1), tile.coords(:,2), 10*ones(size(tile.coords,1),1), '.k', 'markersize', 7)
            set(gca, 'Color', [0.8 0.8 0.8])
            %set(gca, 'xtick', [0.5 1])
            %set(gca, 'ytick', [0 1 2])
            
            % properties
            shading interp; colormap('hot'); colormap(colormap); 
            h = colorbar; 
            %ylabel(h,'$$C_p$$','interpreter','latex')
            %pbaspect([1 2 1])
            %title('$$C_{p,mean}$$','interpreter','latex')
            caxis(axLim); 
            xlabel('$$x[m]$$','interpreter','latex')
            axis([0 1 0.8 2])
            ylabel('$$y[m]$$','interpreter','latex')
            set(gca, 'fontsize', 18)  
            %axis([0.015 1 0.8 2])
            %xlim([0 1]); ylim([0.8462 2])
            
            fig1=gca;
            
% optional outputs
varargout{1} = axLim;
varargout{2} = fig1;

