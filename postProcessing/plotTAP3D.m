function plotTAP3D(tile, taps, angle)

hold on
plotcube([1 0.3 2], [0 0 0], .8, [0.8 0.8 0.8])
pbaspect([1 0.5 2])
set(gca,'xtick',[]); set(gca,'xticklabel',[])
set(gca,'ytick',[]); set(gca,'yticklabel',[])
set(gca,'ztick',[]); set(gca,'zticklabel',[])
if angle == 0 || angle == 180
    a = annotation('textarrow', [0.25 0.35], [0.4 0.45], 'String','U ');
    a.FontSize = 20;
elseif angle == 20
    a = annotation('textarrow', [0.25 0.35], [0.5 0.5], 'String','U ');
    a.FontSize = 20;    
end

if length(taps) == 1    

    plot3(tile.coords(:,1), -0.01*ones(length(tile.coords(:,1)),1), tile.coords(:,2), '.r', 'markersize', 20)
    ylim([-0.1 0.3])
    axis off
   
else
    [~, ~, index] = intersect(taps, tile.taps(:,1:5), 'rows', 'stable');
    
    plot3(tile.coords(:,1), zeros(size(tile.taps,1),1), tile.coords(:,2), '.k', 'markersize', 7)
    plot3(tile.coords(index,1), -0.01*ones(length(index),1), tile.coords(index,2), '.r', 'markersize', 20)
    ylim([-0.1 0.3])
    axis off
end