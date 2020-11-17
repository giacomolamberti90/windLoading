function pressureInstantaneousContour(tile, timeIndex, cLim, list)

    % Function that reads time step data structure of tile A/B or D 
    % the axis limits (optional)
    % and returns: time history
    
% switch nargin
%     case 2
%         % axis limit
%         axLim(1) = min(tile.timeHistory(t,:));
%         axLim(2) = max(tile.timeHistory(t,:)); 
% end

% find pressure peak
[~, ~, index] = intersect(list, tile.taps(:,1:5), 'rows', 'stable');

% Triangulation
DT = delaunayTriangulation(double(tile.coords));

% variable
var = tile.timeHistory(timeIndex,:)';

        hold on
            % contour plot
            patch('Faces', DT.ConnectivityList, 'Vertices', [tile.coords, var],...
                  'FaceVertexCData', var, 'EdgeColor', 'none')
            
            % plot taps
            plot3(tile.coords(:,1), tile.coords(:,2), ones(size(tile.coords,1),1), '.k', 'markersize', 7)
            
            % properties
            title('$$C_p$$','interpreter','latex')            
            shading interp; colormap('hot'); colormap(colormap); colorbar; caxis(cLim)
            %pbaspect([1 1 1])
            xlabel('$$x[m]$$','interpreter','latex'); ylabel('$$y[m]$$','interpreter','latex'); %caxis(cLim)
            set(gca, 'fontsize', 20)
            
%             for i = 1:size(list,1)
%                 plot3(tile.coords(index(i),1), tile.coords(index(i),2), 10, 'o', 'linewidth', 3)
%             end
