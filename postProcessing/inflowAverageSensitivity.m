function [tile_min, tile_mean, tile_max] = inflowAverageSensitivity(tile, stat, param)

    % Function that read LES data structures an returns:
        
        % structure containing:
        %         taps: name and order of pressure taps            [Ntaps,1]        
        %       coords: coordinates of pressure taps               [Ntaps,2]
        %         mean: mean of pressure coefficients              [Ntaps,1]
        %          std: root mean square of pressure coefficients  [Ntaps,1]
        %    CI95_mean: confidence bounds on the mean              [Ntaps,2]
        %     CI95_std: confidence bounds on the standard dev      [Ntaps,2]
        
%% LES
% sensitivity analysis
for i = 1:length(tile)
    
    if stat == 'mean'
        variable(:,i) = tile{i}.mean;
    end
    
    if stat == 'stan'
        variable(:,i) = tile{i}.std;
    end
    
    if stat == 'peak'
        variable(:,i) = tile{i}.peak;
    end
    
    if stat == 'desi'
        [variable(:,i), area] = designPressure('no', tile{i}, 180, 6, 0.22, 'cook');
    end
end

%% "Total Sobol indices"
allIndex = [];
if param == 'z_0'
    %z0
    S_min = []; S_mean = []; S_max = [];
    for i = 1:9
        index  = [3*(i-1) + 1, 3*(i-1) + 2, 3*(i-1) + 3]
        S_min  = [S_min,  variable(:,index(1))];
        S_mean = [S_mean, variable(:,index(2))];
        S_max  = [S_max,  variable(:,index(3))];
        allIndex = [allIndex; index];
    end
elseif  param == 'tke'
    %k
    S_min = []; S_mean = []; S_max = [];
    for i = 1:3
        for j = 1:3
            index = [i + (j-1)*9, i+3 + (j-1)*9, i+6 + (j-1)*9]
            S_min  = [S_min,  variable(:,index(1))];
            S_mean = [S_mean, variable(:,index(2))];
            S_max  = [S_max,  variable(:,index(3))];
            allIndex = [allIndex; index];            
        end
    end
elseif param == 'T_u'
    %Tu
    S_min = []; S_mean = []; S_max = [];
    for i = 1:9
        index = [i, i+9, i+18]
        S_min  = [S_min,  variable(:,index(1))];
        S_mean = [S_mean, variable(:,index(2))];
        S_max  = [S_max,  variable(:,index(3))];
        allIndex = [allIndex; index];
    end
end

%% Structures
tile_min  = tile{14};
tile_mean = tile{14};
tile_max  = tile{14};

% figure
% for k = 1:9
%     
%     subplot(3,3,k)
%     
%     tile_min.peak  = S_min(:,k);
%     tile_mean.peak = S_mean(:,k);
%     tile_max.peak  = S_max(:,k);
%     
%     pressureProfile_all(tile_min,  'A', stat, '-ob');
%     pressureProfile_all(tile_mean, 'A', stat, '-ok');
%     pressureProfile_all(tile_max,  'A', stat, '-or');
% 
% %     semilogx(area, S_min(:,k), '-ob')
% %     hold on
% %     semilogx(area, S_mean(:,k), '-ok')
% %     semilogx(area, S_max(:,k), '-or')
% %     set(gca,'fontsize',22)
% %     set(gca, 'XScale', 'log');
%     
%     tit = sprintf('%d, %d, %d', allIndex(k,:));
%     title(tit);
%     
%     if param == 'z_0'
%         legend({'$$z_0$$ - min','$$z_0$$ - mean','$$z_0$$ - max'}, 'interpreter', 'latex', 'Location', 'NorthWest')
%     elseif param  == 'tke'
%         legend({'$$k$$ - min','$$k$$ - mean','$$k$$ - max'}, 'interpreter', 'latex', 'Location', 'NorthWest')
%     elseif param  == 'T_u'
%         legend({'$$T_u$$ - min','$$T_u$$ - mean','$$T_u$$ - max'}, 'interpreter', 'latex', 'Location', 'NorthWest')
%     end
% end

if stat == 'mean'
    tile_min.mean  = mean(S_min');
    tile_mean.mean = mean(S_mean');
    tile_max.mean  = mean(S_max');
end

if stat == 'stan'
    tile_min.std  = mean(S_min');
    tile_mean.std = mean(S_mean');
    tile_max.std  = mean(S_max');
end

if stat == 'peak'
    tile_min.peak  = mean(S_min');
    tile_mean.peak = mean(S_mean');
    tile_max.peak  = mean(S_max');
end

if stat == 'desi'
    tile_min.area    = area;
    tile_min.design  = mean(S_min');
    tile_mean.area   = area;
    tile_mean.design = mean(S_mean');
    tile_max.area    = area;
    tile_max.design  = mean(S_max');
end
