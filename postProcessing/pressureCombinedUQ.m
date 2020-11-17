function [tileA_comb, tileB_comb] = pressureCombinedUQ(tileA_base, tileA_c1c, tileA_c2c, tileA_c3c, ...
                                                       tileB_base, tileB_c1c, tileB_c2c, tileB_c3c)                                                       

    % Function that read .dat file from WoW experiment of tiles A-B 
    % and returns:     
    
        %structure containing:
        %      coords: coordinates of pressure taps               [Ntaps,2]        
        %       names: name and order of pressure taps            [Ntaps,1]
        % timeHistory: time series of pressure coefficients       [Nsamp,Ntaps]
        %        mean: mean of pressure coefficients              [Ntaps,1]
        %         std: root mean square of pressure coefficients  [Ntaps,1]
        
%% define structures
tileA_comb.coords = tileA_base.coords;
tileB_comb.coords = tileB_base.coords;

tileA_comb.taps = tileA_base.taps;
tileB_comb.taps = tileB_base.taps;

tileA_comb.mean = tileA_base.mean;
tileB_comb.mean = tileB_base.mean;

tileA_comb.CI = 'on';
tileB_comb.CI = 'on';

tileA_comb.CI95_mean(:,1) = min([tileA_base.CI95_mean(:,1)'; tileA_c1c.CI95_mean(:,1)'; ...
                            tileA_c2c.CI95_mean(:,1)';  tileA_c3c.CI95_mean(:,1)'])';
                        
tileA_comb.CI95_mean(:,2) = max([tileA_base.CI95_mean(:,2)'; tileA_c1c.CI95_mean(:,2)'; ...
                            tileA_c2c.CI95_mean(:,2)';  tileA_c3c.CI95_mean(:,2)'])';

tileB_comb.CI95_mean(:,1) = min([tileB_base.CI95_mean(:,1)'; tileB_c1c.CI95_mean(:,1)'; ...
                            tileB_c2c.CI95_mean(:,1)';  tileB_c3c.CI95_mean(:,1)'])';
                        
tileB_comb.CI95_mean(:,2) = max([tileB_base.CI95_mean(:,2)'; tileB_c1c.CI95_mean(:,2)'; ...
                            tileB_c2c.CI95_mean(:,2)';  tileB_c3c.CI95_mean(:,2)'])';         
