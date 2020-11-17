function tile = designPressurePanels(tile, dx, dy, window, npeaks, start)

    % Function that reads tile and POE 
    % returns cp design relative to POE, as a function of the area of the
    % panels
        
%%
if min(tile.coords(:,1)) < 0.4;
    
    X_A = 0.0 : dx : 0.6;
    Y_A = 1.4 : dy : 2;
    
    X_B = 0.0  : dx : 0.25;
    Y_B = 0.85 : dy : 1.15;

else
    X_A = 0.4 : dx : 1;
    Y_A = 1.4 : dy : 2;

    X_B = 0.7 : dx : 1;
    Y_B = 0.8 : dy : 1.2;

end


if tile.taps(1) == 'A'
    
    k = 0;
    for j = 1:length(Y_A)-1
       for i = 1:length(X_A)-1
           k = k + 1;
           panels(k,:) = [X_A(i), X_A(i+1), Y_A(j), Y_A(j+1)];
       end
    end
    
elseif tile.taps(1) == 'B'
    
    k = 0;
    for j = 1:length(Y_B)-1
       for i = 1:length(X_B)-1
           k = k + 1;
           panels(k,:) = [X_B(i), X_B(i+1), Y_B(j), Y_B(j+1)];
       end
    end
    
elseif tile.taps(1) == 'D'
    
    k = 0;
    for j = 1:length(Y_D)-1
       for i = 1:length(X_D)-1
           k = k + 1;
           panels(k,:) = [X_D(i), X_D(i+1), Y_D(j), Y_D(j+1)];
       end
    end   
    
end

for k = 1:size(panels,1)
    
    panel = panels(k,:);
    
    % compute area-averaged pressure
    tile.areaAverage = areaAveragedPressure(tile, panel);
    
    if isnan(tile.areaAverage) == 0
        cp_design(k) = gumbel(tile, window, 0.22, npeaks, start);
    else
        cp_design(k) = nan;
    end
end
index = find(isnan(cp_design) == 0);

%cp_design = cp_design(index);
%panels = panels(index,:);

hold on
for i = 1:length(cp_design)
   
    panel = panels(i,:);
    
    X_vertex = [panel(1), panel(2), panel(2), panel(1)];
    Y_vertex = [panel(3), panel(3), panel(4), panel(4)];
    
    if isnan(cp_design(i)) == 1
        patch(X_vertex, Y_vertex, cp_design(i), 'LineStyle', 'none')
    else
        patch(X_vertex, Y_vertex, cp_design(i))
    end
end
colormap(hot)
colorbar

scatter(tile.coords(:,1), tile.coords(:,2), 5, 'k', 'FaceColor', 'flat');

xlabel('$$x[m]$$','interpreter','latex')
ylabel('$$y[m]$$','interpreter','latex')
pbaspect([1 2 1]); %axis([0 1 0 2])
set(gca, 'Color', [0.8 0.8 0.8])
set(gca, 'fontsize', 20)
axis([min(X_A), max(X_A), min(Y_B), max(Y_A)])
caxis([-2.5, -0.5])

box on

tile.design = cp_design;
tile.panels = panels;

