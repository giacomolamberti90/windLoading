function [tileA, tileB] = pressureCharlesLES(path_LES, ang)

    % Function that read LES mean, variance and time series from probes
    % of tiles A-B and returns:
        
        % structure containing:
        %         taps: name and order of pressure taps            [Ntaps,1]        
        %       coords: coordinates of pressure taps               [Ntaps,2]
        %  timeHistory: time series of pressure coefficients       [Ntaps,Nsamp]       
        %         mean: mean of pressure coefficients              [Ntaps,1]
        %          std: root mean square of pressure coefficients  [Ntaps,1]
        
%% taps list - reference
tapsList   = load('/home/giacomol/Desktop/Research/windLoading/windTunnel/matlabCode/ArupGUI3.0/Orifizi-v2.mat');
taps_A_ref = cell2mat(tapsList.nomi{1}');
taps_B_ref = cell2mat(tapsList.nomi{2}');

%% define structures
data_A = dlmread(fullfile(path_LES, sprintf('%i/resultA%i.txt', ang, ang)));
data_B = dlmread(fullfile(path_LES, sprintf('%i/resultB%i.txt', ang, ang)));

% tileA
tileA.CI = 'no';
tileA.taps = taps_A_ref;
tileA.coords = data_A(:,1:2);
tileA.mean = data_A(:,4);
tileA.std = data_A(:,5);
tileA.peak = data_A(:,6)';


% tileB
tileB.CI = 'no';
tileB.taps = taps_B_ref;
tileB.coords = data_B(:,1:2);
tileB.mean = data_B(:,4);
tileB.std = data_B(:,5);
tileB.peak = data_B(:,6)';