function [tileA, tileB, tileD] = pressureLES(path_LES, ang, deltaT, deltaT_filter)

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
taps_D_ref = cell2mat(tapsList.nomi{4}');

path = '/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/';

% taps coordinates
coords_A = load(fullfile(path, 'coords_A0'));
coords_B = load(fullfile(path, 'coords_B0'));
coords_D = tapsList.prese{4}/1000; coords_D(:,2) = 2-coords_D(:,2);

coords_A_ref = coords_A;
coords_B_ref = coords_B;

if ang == 180
    coords_A = load(fullfile(path, 'coords_A180'));
    coords_B = load(fullfile(path, 'coords_B180'));

    coords_A_ref = coords_A;
    coords_B_ref = coords_B;

elseif ang == 20
    coords_U = load(fullfile(path, 'coords_20_upper'));
    coords_A = load(fullfile(path, 'coords_A20'));
    coords_B = load(fullfile(path, 'coords_B20'));
end

%% data
% read probesls
Nprobes_V = 7; NV = 1;

data_A = [];
data_B = [];

if ang == 0
    Nprobes_A = 224;
    Nprobes_B = 223;
    % area-averaged panels
    panelA = [0.96 1.00 1.94 2.00];
    panelB = [0.96 1.00 0.97 1.03];
    dire = dir(fullfile(path_LES, sprintf('/postProcessing/A%d', ang)));
elseif ang == 180
    Nprobes_A = 224;
    Nprobes_B = 223;
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];
    dire = dir(fullfile(path_LES, sprintf('/postProcessing/A%d', ang)));
else
    Nprobes_U = size(coords_U, 1);
    Nprobes_A = 224;
    Nprobes_B = 223;    
    % area-averaged panels
    panelA = [0.96 1.00 1.94 2.00];
    dire = dir(fullfile(path_LES, sprintf('/postProcessing/upper_%d', ang)));
end

% last time-step
dire_last = dir(fullfile(path_LES, '/postProcessing/probes_x_2m_old'));
for i = 1:length(dire_last)-2
    tlast(i) = str2num(dire_last(2+i).name);
end
tlast = sort(tlast);
T_end = tlast(end);

% time in probes
for i = 1:length(dire)-2
    tprobes(i) = str2num(dire(2+i).name);
end
tprobes = sort(tprobes);
tprobes = tprobes(tprobes <= T_end);

% N = 1;
N = length(tprobes)-1;

for i = length(tprobes)-N:length(tprobes)
    
    tprobes(i)

    if ang == 0 || ang == 180
        
        Nprobes = Nprobes_A;
        
        % time series
        fid    = fopen(fullfile(path_LES, sprintf('postProcessing/A%d/%s/p', ang, num2str(tprobes(i)))));
        format = repmat('%f', 1, Nprobes+1);
        data   = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes+2);
        fclose(fid); 

        data_A = [data_A; cell2mat(data)];
        clear data
    
        fid    = fopen(fullfile(path_LES, sprintf('postProcessing/B%d/%s/p', ang, num2str(tprobes(i)))));
        format = repmat('%f', 1, Nprobes_B+1);
        data   = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_B+2);
        fclose(fid);

        data_B = [data_B; cell2mat(data)];
        clear data
    else
        
        Nprobes = Nprobes_U;
        
        % time series
        fid    = fopen(fullfile(path_LES, sprintf('postProcessing/upper_%d/%s/p', ang, num2str(tprobes(i)))));
        format = repmat('%f', 1, Nprobes+1);
        data   = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes+2);
        fclose(fid); 

        data_A = [data_A; cell2mat(data)];
        clear data        
    end
end

% fid    = fopen(fullfile(path_LES, sprintf('postProcessing/probes_x/%s/UMean', num2str(tprobes(end)))));
% format = repmat('%f', 1, 3*(Nprobes_V+1));
% data_U = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_V+2);
% fclose(fid);
% 
% U_1m = data_U{2+3*(NV-1)}(end);
% sprintf('Uref = %f', U_1m)

fid    = fopen(fullfile(path_LES, sprintf('postProcessing/probes_x_2m_old/%s/UMean', num2str(tprobes(end)))));
format = repmat('%f', 1, 3*(Nprobes_V+1));
data_U = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_V+2);
fclose(fid);

U_2m = data_U{2+3*(NV-1)}(end);
U = U_2m;

% pressure time series
[time_A, it_A, ~] = unique(data_A(:,1));
    
if ang == 0 || ang == 180
    
    [time_B, it_B, ~] = unique(data_B(:,1));
    
    for i = 1:Nprobes_A
        p_A(:,i) = data_A(it_A,i+1)/(0.5*U^2);
        if p_A(:,i) <= -1e6
            p_A(:,i) = nan;
        end
    end    

    for i = 1:Nprobes_B
        p_B(:,i) = data_B(it_B,i+1)/(0.5*U^2);
        if p_B(:,i) <= -1e6
            p_B(:,i) = nan;
        end
    end        
    
    % mean 
    fid     = fopen(fullfile(path_LES, sprintf('postProcessing/A%d/%s/pMean', ang, num2str(tprobes(end)))));
    format  = repmat('%f', 1, Nprobes_A+1);
    datam_A = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_A+2);
    fclose(fid);

    fid     = fopen(fullfile(path_LES, sprintf('postProcessing/B%d/%s/pMean', ang, num2str(tprobes(end)))));
    format  = repmat('%f', 1, Nprobes_B+1);
    datam_B = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_B+2);
    fclose(fid);

    %fid     = fopen(fullfile(path_LES, sprintf('postProcessing/D0/%s/pMean', ang, num2str(tprobes(end)))));
    %format  = repmat('%f', 1, Nprobes_B+1);
    %datam_D = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_B+2);
    %fclose(fid);
    
    % variance 
    fid     = fopen(fullfile(path_LES, sprintf('postProcessing/A%d/%s/pPrime2Mean', ang, num2str(tprobes(end)))));
    format  = repmat('%f', 1, Nprobes_A+1);
    datav_A = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_A+2);
    fclose(fid); 

    fid     = fopen(fullfile(path_LES, sprintf('postProcessing/B%d/%s/pPrime2Mean', ang, num2str(tprobes(end)))));
    format  = repmat('%f', 1, Nprobes_B+1);
    datav_B = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_B+2);
    fclose(fid);

    %fid     = fopen(fullfile(path_LES, sprintf('postProcessing/D%d/%s/pPrime2Mean', ang, num2str(tprobes(end)))));
    %format  = repmat('%f', 1, Nprobes_B+1);
    %datav_D = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_B+2);
    %fclose(fid);
    
    % Pressure mean and root mean square
    % tile A
    for i = 1:Nprobes_A
        cp_mean_A(i) = datam_A{i+1}(end)/(0.5*U^2);
        cp_std_A(i)  = sqrt(datav_A{i+1}(end))/(0.5*U^2);
    end

    % tile B-D
    for i = 1:Nprobes_B
        cp_mean_B(i) = datam_B{i+1}(end)/(0.5*U^2);
        cp_std_B(i)  = sqrt(datav_B{i+1}(end))/(0.5*U^2);
        
        %cp_mean_D(i) = datam_D{i+1}(end)/(0.5*U^2);
        %cp_std_D(i)  = sqrt(datav_D{i+1}(end))/(0.5*U^2);        
    end   
    
else
    idx  = [];
    for i = 1:Nprobes_U
        p_A(:,i) = data_A(it_A,i+1)/(0.5*U^2);
        if p_A(:,i) <= -1e6
            p_A(:,i) = nan;
        else
            idx = [idx; i];            
        end
    end    
    
    % mean
    cp_mean_A = dlmread(fullfile(path_LES, sprintf('cp_mean_A%d.out', ang)))';
    cp_mean_B = dlmread(fullfile(path_LES, sprintf('cp_mean_B%d.out', ang)))';
%     
    % std
    cp_std_A = dlmread(fullfile(path_LES, sprintf('cp_rms_A%d.out', ang)))';
    cp_std_B = dlmread(fullfile(path_LES, sprintf('cp_rms_B%d.out', ang)))';

%     % mean 
%     fid     = fopen(fullfile(path_LES, sprintf('postProcessing/upper_%d/%s/pMean', ang, num2str(tprobes(end)))));
%     format  = repmat('%f', 1, Nprobes_U+1);
%     datam_A = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_U+2);
%     fclose(fid);
%     
%     % variance 
%     fid     = fopen(fullfile(path_LES, sprintf('postProcessing/upper_%d/%s/pPrime2Mean', ang, num2str(tprobes(end)))));
%     format  = repmat('%f', 1, Nprobes_U+1);
%     datav_A = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_U+2);
%     fclose(fid);
%     
%     % tile A
%     for i = 1:Nprobes_U
%         if datam_A{i+1}(end)
%             cp_mean_A(i) = datam_A{i+1}(end)/(0.5*U^2);
%             cp_std_A(i)  = sqrt(datav_A{i+1}(end))/(0.5*U^2);
%         end
%     end
end

%% define structures
step = floor(deltaT_filter/deltaT);

% tileA
tileA.CI = 'no';
tileA.taps = taps_A_ref;
if ang == 20
    tileA.coords = coords_U(idx,:);
else
    tileA.coords = coords_A_ref;
end
% tileA.coords = coords_A_ref;
tileA.mean = cp_mean_A';
tileA.std = cp_std_A';
tileA.time = time_A(1:step:end)-time_A(1);
tileA.timeHistory = p_A(1:step:end, :);
% tileA.U_1m = U_1m;
tileA.U_2m = U_2m;
tileA.areaAverage = areaAveragedPressure(tileA, panelA);
% tileA.peak = pressurePeak(tileA, 6, 0.22, 'cook');

%tileB
tileB = {};
tileB.CI = 'no';
tileB.taps = taps_B_ref;
tileB.coords = coords_B_ref;
tileB.mean = cp_mean_B';
tileB.std = cp_std_B';
% tileB.U_1m = U_1m;
tileB.U_2m = U_2m;
if ang == 0 || ang == 180
    tileB.time = time_B(1:step:end)-time_B(1);
    tileB.timeHistory = p_B(1:step:end,:);    
    tileB.areaAverage = areaAveragedPressure(tileB, panelB);
%     tileB.peak = pressurePeak(tileB, 6, 0.22, 'cook');
end

%tileD
tileD.CI = 'no';
tileD.taps = taps_D_ref;
tileD.coords = coords_D;
% if ang == 0 || ang == 180
%     tileD.mean = pmean_D'/qD;
%     tileD.std  = pstd_D'/qD;
%     tileD.U    = U;
% end
