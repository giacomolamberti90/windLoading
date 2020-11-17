function [tileA_sens, tileB_sens, tileD_sens] = pressureCombinedLES(tileA, tileB, tileD)

    % Function that read LES data structures an returns:
        
        % structure containing:
        %         taps: name and order of pressure taps            [Ntaps,1]        
        %       coords: coordinates of pressure taps               [Ntaps,2]
        %         mean: mean of pressure coefficients              [Ntaps,1]
        %          std: root mean square of pressure coefficients  [Ntaps,1]
        %    CI95_mean: confidence bounds on the mean              [Ntaps,2]
        %     CI95_std: confidence bounds on the standard dev      [Ntaps,2]
        
%% LES
% for i = 1:27
%     sprintf('workdir = %i', i)
%     [tileA{i}, tileB{i}] = pressureLES(fullfile(path_SENS, sprintf('workdir.%i',i)), ang, 15);
% end
    
% sensitivity analysis
for i = 1:length(tileA)
    
    % extract mean
    meanA(:,i) = tileA{i}.mean;
    meanB(:,i) = tileB{i}.mean;
    %meanD(:,i) = tileD{i}.mean;
    
    % extract rms
    rmsA(:,i) = tileA{i}.std;
    rmsB(:,i) = tileB{i}.std;
    %stdD(:,i) = tileD{i}.std;
    
    % extract peak
    peakA(:,i) = tileA{i}.peak;
    peakB(:,i) = tileB{i}.peak;
end

%% "Total Sobol indices" 
%z0
for i = 1:9
    
    index = [3*(i-1) + 1, 3*(i-1) + 2, 3*(i-1) + 3];
    S1_mean_A(:,i) = var(meanA(:,index)')./var(meanA');
    S1_mean_B(:,i) = var(meanB(:,index)')./var(meanB');      
    S1_rms_A(:,i) = var(rmsA(:,index)')./var(rmsA');
    S1_rms_B(:,i) = var(rmsB(:,index)')./var(rmsB');    
    S1_peak_A(:,i) = var(peakA(:,index)')./var(peakA');
    S1_peak_B(:,i) = var(peakB(:,index)')./var(peakB');
end
S1_mean_A = mean(S1_mean_A');
S1_mean_B = mean(S1_mean_B');
S1_rms_A = mean(S1_rms_A');
S1_rms_B = mean(S1_rms_B');
S1_peak_A = mean(S1_peak_A');
S1_peak_B = mean(S1_peak_B');

%k
k = 1;
for i = 1:3
    for j = 1:3
        index = [i + (j-1)*9, i+3 + (j-1)*9, i+6 + (j-1)*9];
        S2_mean_A(:,i) = var(meanA(:,index)')./var(meanA');
        S2_mean_B(:,i) = var(meanB(:,index)')./var(meanB');          
        S2_rms_A(:,k) = var(rmsA(:,index)')./var(rmsA');
        S2_rms_B(:,k) = var(rmsB(:,index)')./var(rmsB');  
        S2_peak_A(:,k) = var(peakA(:,index)')./var(peakA');
        S2_peak_B(:,k) = var(peakB(:,index)')./var(peakB');
        k = k+1;
    end
end
S2_mean_A = mean(S2_mean_A');
S2_mean_B = mean(S2_mean_B');
S2_rms_A = mean(S2_rms_A');
S2_rms_B = mean(S2_rms_B');
S2_peak_A = mean(S2_peak_A');
S2_peak_B = mean(S2_peak_B');

%Tu
for i = 1:9
    index = [i, i+9, i+18];
    S3_mean_A(:,i) = var(meanA(:,index)')./var(meanA');
    S3_mean_B(:,i) = var(meanB(:,index)')./var(meanB');      
    S3_rms_A(:,i) = var(rmsA(:,index)')./var(rmsA');
    S3_rms_B(:,i) = var(rmsB(:,index)')./var(rmsB');  
    S3_peak_A(:,i) = var(peakA(:,index)')./var(peakA');
    S3_peak_B(:,i) = var(peakB(:,index)')./var(peakB');
end
S3_mean_A = mean(S3_mean_A');
S3_mean_B = mean(S3_mean_B');
S3_rms_A = mean(S3_rms_A');
S3_rms_B = mean(S3_rms_B');
S3_peak_A = mean(S3_peak_A');
S3_peak_B = mean(S3_peak_B');


%% "Sobol indices" 
% %z0
% for i = 1:3    
%     index = i:3:27;
%     E1_mean_A(:,i) = mean(meanA(:,index)');
%     E1_mean_B(:,i) = mean(meanB(:,index)');        
%     E1_rms_A(:,i) = mean(rmsA(:,index)');
%     E1_rms_B(:,i) = mean(rmsB(:,index)');    
%     E1_peak_A(:,i) = mean(peakA(:,index)');
%     E1_peak_B(:,i) = mean(peakB(:,index)');
% end
% S1_mean_A = var(E1_mean_A')./var(meanA');
% S1_mean_B = var(E1_mean_B')./var(meanB');
% S1_rms_A = var(E1_rms_A')./var(rmsA');
% S1_rms_B = var(E1_rms_B')./var(rmsB');
% S1_peak_A = var(E1_peak_A')./var(peakA');
% S1_peak_B = var(E1_peak_B')./var(peakB');
% 
% %k
% k = 1;
% for i = 1:3:7
%     index = [i:i+2, i+9:i+11, i+18:i+20];
%     E2_mean_A(:,k) = mean(meanA(:,index)');
%     E2_mean_B(:,k) = mean(meanB(:,index)');            
%     E2_rms_A(:,k) = mean(rmsA(:,index)');
%     E2_rms_B(:,k) = mean(rmsB(:,index)');  
%     E2_peak_A(:,k) = mean(peakA(:,index)');
%     E2_peak_B(:,k) = mean(peakB(:,index)');
%     k = k+1;
% end
% S2_mean_A = var(E2_mean_A')./var(meanA');
% S2_mean_B = var(E2_mean_B')./var(meanB');
% S2_rms_A = var(E2_rms_A')./var(rmsA');
% S2_rms_B = var(E2_rms_B')./var(rmsB');
% S2_peak_A = var(E2_peak_A')./var(peakA');
% S2_peak_B = var(E2_peak_B')./var(peakB');
% 
% %Tu
% k = 1;
% for i = 1:9:19
%     index = i:i+8;
%     E3_mean_A(:,k) = mean(meanA(:,index)');
%     E3_mean_B(:,k) = mean(meanB(:,index)');            
%     E3_rms_A(:,k) = mean(rmsA(:,index)');
%     E3_rms_B(:,k) = mean(rmsB(:,index)');  
%     E3_peak_A(:,k) = mean(peakA(:,index)');
%     E3_peak_B(:,k) = mean(peakB(:,index)');
%     k = k+1;
% end
% S3_mean_A = var(E3_mean_A')./var(meanA');
% S3_mean_B = var(E3_mean_B')./var(meanB');
% S3_rms_A = var(E3_rms_A')./var(rmsA');
% S3_rms_B = var(E3_rms_B')./var(rmsB');
% S3_peak_A = var(E3_peak_A')./var(peakA');
% S3_peak_B = var(E3_peak_B')./var(peakB');

%% define structures
% tile A
tileA_sens.CI          = 'on';
tileA_sens.taps        = tileA{1}.taps;
tileA_sens.coords      = tileA{1}.coords;
tileA_sens.mean        = meanA(:,14);
tileA_sens.std         = rmsA(:,14);
tileA_sens.CI95_mean   = [min(meanA')', max(meanA')'];
tileA_sens.CI95_std    = [min(rmsA')',  max(rmsA')'];
tileA_sens.CI95_peak   = [min(peakA')', max(peakA')'];
% tileA_sens.sobol       = [S1_mean_A', S2_mean_A', S3_mean_A'];
tileA_sens.sobol       = [S1_rms_A', S2_rms_A', S3_rms_A'];
% tileA_sens.sobol       = [S1_peak_A', S2_peak_A', S3_peak_A'];
tileA_sens.peak        = peakA(:,14)';

% tile B
tileB_sens.CI          = 'on';
tileB_sens.taps        = tileB{1}.taps;
tileB_sens.coords      = tileB{1}.coords;
tileB_sens.mean        = meanB(:,14);
tileB_sens.std         = rmsB(:,14);
tileB_sens.CI95_mean   = [min(meanB')', max(meanB')'];
tileB_sens.CI95_std    = [min(rmsB')',  max(rmsB')'];
tileB_sens.CI95_peak   = [min(peakB')', max(peakB')'];
% tileB_sens.sobol       = [S1_mean_B', S2_mean_B', S3_mean_B'];
tileB_sens.sobol       = [S1_rms_B', S2_rms_B', S3_rms_B'];
% tileB_sens.sobol       = [S1_peak_B', S2_peak_B', S3_peak_B'];
tileB_sens.peak        = peakB(:,14)';

% tile D
tileD_sens.CI          = 'on';
tileD_sens.taps        = tileD{1}.taps;
tileD_sens.coords      = tileD{1}.coords;
%tileD_sens.mean        = mean(meanD,2);
%tileD_sens.std         = mean(stdD,2);
%tileD_sens.CI95_mean   = [min(meanD')', max(meanD')'];
% tileD_sens.CI95_std    = [min(stdD')', max(stdD')'];
