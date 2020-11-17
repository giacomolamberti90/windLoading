function [tileA, tileB, tileD] = pressureWoW_HR(angle)

    % Function that read .dat file from WoW experiment of tiles A-B-D 
    % and returns:     
    
        %structure containing:
        %      coords: coordinates of pressure taps               [Ntaps,2]        
        %       names: name and order of pressure taps            [Ntaps,1]
        % timeHistory: time series of pressure coefficients       [Nsamp,Ntaps]
        %        mean: mean of pressure coefficients              [Ntaps,1]
        %         std: root mean square of pressure coefficients  [Ntaps,1]
        
%% taps list - reference
tapsList   = load('/home/giacomol/Desktop/Research/windLoading/windTunnel/matlabCode/ArupGUI3.0/Orifizi-v2.mat');
coords_A   = tapsList.prese{1};
coords_B   = tapsList.prese{2};
coords_D   = tapsList.prese{4};
taps_A_ref = cell2mat(tapsList.nomi{1}');
taps_B_ref = cell2mat(tapsList.nomi{2}');
taps_D_ref = taps_B_ref;

%% data
fsamp   = 520;
Ntaps_A = 224;
Ntaps_B = 223; % last one 'not used'

% dynamic pressure
path_q = '/home/storage/WoW_data/Raw_Data_Untouched/Cobra_Probe_Data/Braced';
file_q = '2017-102n_O_A_AL_HR1_CB_1to50_055perc_000deg_Test20 (Ve).thA';
[u, ~, ~, ~, Settings] = ReadTHFile(fullfile(path_q,file_q));

% reference quantities
U = mean(u);
T = Settings.Tmean + 273.17;
P = Settings.Pbaro;
R = 287;
q = 0.5*1.16*U^2;%0.5*(P/(R*T))*U^2;

path = '/home/storage/WoW_data/Raw_Data_Untouched/Scanivalve_Data/Braced';

if angle == 0 %(WoW 180)

elseif angle == 10 %(WoW 170)
    fileA_HR1 = '2017-102n_O_A_AL_HR1_1to50_055perc_170deg_Test21.dat';
    fileA_HR2 = '2017-102n_O_A_AL_HR2_1to50_055perc_170deg_Test22.dat';
elseif angle == 20 %(WoW 160)
    fileA_HR1 = '2017-102n_O_A_AL_HR1_1to50_055perc_160deg_Test21.dat';
    fileA_HR2 = '2017-102n_O_A_AL_HR2_1to50_055perc_160deg_Test22.dat';
elseif angle == 180 %(WoW 0)

end
dataA_HR1 = fullfile(path,fileA_HR1);
dataA_HR2 = fullfile(path,fileA_HR2);

% pressure data
psi2Pa  = 6894.76;
pdata_A_HR1 = psi2Pa*importScanivalveDataFile(dataA_HR1);
pdata_A_HR2 = psi2Pa*importScanivalveDataFile(dataA_HR2);

pdata_A = [pdata_A_HR1, pdata_A_HR2];

% % transfer function
% path_TF   = '/home/giacomol/Desktop/Research/windLoading/windTunnel/WoW/transferFunction';
% file_TF_A = '2017-102n-Stanford-High-rise_FS_1652Hz_TS13_5in_plus_6in_CS12in_TFunc_Test11_15min_FF_250Hz_FSO_520Hz.mat';
% file_TF_B = '2017-102n-Stanford-High-rise_FS_1652Hz_TS77_25in_plus_6in_CS12in_TFunc_Test12_15min_FF_250Hz_FSO_520Hz.mat';
% TFFName_A = fullfile(path_TF,file_TF_A);
% TFFName_B = fullfile(path_TF,file_TF_B);

% % apply transfer function
% pdata_A = applyTF(pdata_A, TFFName_A, 520, 250);
% pdata_B = applyTF(pdata_B, TFFName_B, 520, 250);
% pdata_D = applyTF(pdata_D, TFFName_B, 520, 250);

% mean pressure from Z taps
pmeanZ_A = 0;%mean(nanmean(pdata_A(:,Ntaps_A+1:end)));
pmeanZ_B = 0;%mean(nanmean(pdata_B(:,Ntaps_B+1:end)));
pmeanZ_D = 0;%mean(nanmean(pdata_D(:,Ntaps_B+1:end)));

% taps names
[~,txt,~] = xlsread('/home/storage/WoW_data/TestDocuments/Drawings/2017-102n-Tubing-Table.xlsx', 'Tile A');
taps_A    = txt(2:end,2);

% remove corrupted data 
%[~, ~, index] = setxor({'A0710'}, taps_A(1:Ntaps_A,:));
[~, ~, index] = setxor({''}, taps_A(1:Ntaps_A,:));

taps_A     = taps_A(index,:);
taps_A_ref = taps_A_ref(index,:);
coords_A   = coords_A(index,:);
pdata_A    = pdata_A(:,index);

% convert to array
taps_A = cell2mat(taps_A);

% re-arrange data in reference order
[~, index_A_ref, index_A] = intersect(taps_A_ref, taps_A(:,1:5), 'rows', 'stable');

pdata_A = pdata_A(:,index_A_ref);

% reorder coordinates and convert to meters
coords_A = coords_A(index_A_ref,:)/1000;

% convert coordinates to CFD reference system
coords_A(:,1) = 1+coords_A(:,1);

if angle == 180;
    coords_A(:,1) = 1-coords_A(:,1);  
end

%% Pressure mean and variance
% tile A
pmean_A = nanmean(pdata_A,1);
pstd_A  = nanstd(pdata_A,1);

% prssure coefficients
cp_A     = pdata_A/q;
cpmean_A = (pmean_A-pmeanZ_A)/q;
cpstd_A  = pstd_A/q;

% define structures
% tileA
tileA.CI          = 'no';
tileA.taps        = taps_A(index_A_ref,:);
tileA.coords      = coords_A;
tileA.time        = linspace(0, size(pdata_A,1)/fsamp, size(pdata_A,1));
tileA.timeHistory = cp_A;
tileA.mean        = cpmean_A';
tileA.std         = cpstd_A';
% tileB
tileB = [];
tileD = [];
