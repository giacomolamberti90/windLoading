function tile_comb = pressureCombinedUQ_all(tile_base, tile_c1c, tile_c2c, tile_c3c)                                                       

    % Function that read .dat file from WoW experiment of tiles A-B 
    % and returns:     
    
        %structure containing:
        %      coords: coordinates of pressure taps               [Ntaps,2]        
        %       names: name and order of pressure taps            [Ntaps,1]
        % timeHistory: time series of pressure coefficients       [Nsamp,Ntaps]
        %        mean: mean of pressure coefficients              [Ntaps,1]
        %         std: root mean square of pressure coefficients  [Ntaps,1]
        
%% define structures
tile_comb.coords = tile_base.coords;
tile_comb.taps = tile_base.taps;
tile_comb.mean = tile_base.mean;
tile_comb.CI = 'on';

tile_comb.CI95_mean(:,1) = min([tile_base.CI95_mean(:,1)'; tile_c1c.CI95_mean(:,1)'; ...
                                 tile_c2c.CI95_mean(:,1)';  tile_c3c.CI95_mean(:,1)'])';
                        
tile_comb.CI95_mean(:,2) = max([tile_base.CI95_mean(:,2)'; tile_c1c.CI95_mean(:,2)'; ...
                                 tile_c2c.CI95_mean(:,2)';  tile_c3c.CI95_mean(:,2)'])';