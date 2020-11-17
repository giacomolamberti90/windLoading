function plotTAP(tile, taps)

hold on
set(gca, 'fontsize', 20)
pbaspect([1 2 1])
xlabel('x[m]')
ylabel('y[m]')
axis([0 1 0 2])
set(gca,'Color', [0.8 0.8 0.8])

if length(taps) == 1 
    var = 0.75*ones(size(tile.coords,1), 1);
    DT  = delaunayTriangulation(double(tile.coords));

    hold on
        % contour plot                
        patch('Faces', DT.ConnectivityList, 'Vertices', [tile.coords, var],...
              'FaceVertexCData', var, 'EdgeColor', 'none','FaceAlpha', 0.5);
          
   shading interp; colormap('hot'); caxis([0,2]);
   
else
    [~, ~, index] = intersect(taps, tile.taps(:,1:5), 'rows', 'stable');
    
    plot(tile.coords(:,1), tile.coords(:,2), '.k', 'markersize', 7)
    plot(tile.coords(index,1), tile.coords(index,2), '.r', 'markersize', 15)
end