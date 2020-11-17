function [tileA, tileB, tileD] = pressureWoW(angle)

    % Function that read .dat file from WoW experiment of tiles A-B-D 
    % and returns:     
    
        % structure containing:
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
taps_D_ref = cell2mat(tapsList.nomi{4}');

%% data
fsamp   = 520;
Ntaps_A = 224;
Ntaps_B = 223; % last one 'not used'

% dynamic pressure
path_q = '/home/storage/WoW_data/Raw_Data_Untouched/Cobra_Probe_Data/Braced';
file_q = '2017-102n_O_A_BR_CB_1to50_018perc_000deg_Test1 (Ve).thA';
% file_q = '2017-102n_O_B_BR_1to50_018perc_000deg_Test3 (Ve).thB';
[u, ~, ~, ~, Settings] = ReadTHFile(fullfile(path_q,file_q));

z0    = 0.025;
ustar = 0.41*mean(u)/log(2/z0);

% reference quantities
U  = ustar/0.41*log(2/z0); %9.26
T  = Settings.Tmean + 273.17;
P  = Settings.Pbaro;
R  = 287;
qA = 0.5*(P/(R*T))*U^2; qB = qA; qD = qA;

path = '/home/storage/WoW_data/Raw_Data_Untouched/Scanivalve_Data/Braced';

if angle == 0 %(WoW 180)
    fileA = '2017-102n_O_A_BR_1to50_018perc_180deg_Test2.dat';
    fileB = '2017-102n_O_B_BR_1to50_018perc_180deg_Test4.dat'; 
    fileD = '2017-102n_O_D_BR_1to50_018perc_180deg_Test6.dat';
    % area-averaged panels
    panelA = [0.96 1.00 1.94 2.00];
    panelB = [0.96 1.00 0.97 1.03];
elseif angle == 10 %(WoW 170)
    fileA = '2017-102n_O_A_BR_1to50_018perc_170deg_Test2.dat';
    %fileA = '2017-102n_O_A_AL_1to50_018perc_170deg_Test17.dat';
    fileB = '2017-102n_O_B_BR_1to50_018perc_170deg_Test4.dat';  
    fileD = '2017-102n_O_D_BR_1to50_018perc_170deg_Test6.dat';
elseif angle == 20 %(WoW 160)
    fileA = '2017-102n_O_A_BR_1to50_018perc_160deg_Test2.dat';
    %fileA = '2017-102n_O_A_AL_1to50_018perc_160deg_Test17.dat';
    fileB = '2017-102n_O_B_BR_1to50_018perc_160deg_Test4.dat';    
    fileD = '2017-102n_O_D_BR_1to50_018perc_160deg_Test6.dat';
    % area-averaged panels
    panelA = [0.96 1.00 1.94 2.00];
    panelB = [0.96 1.00 0.97 1.03];
elseif angle == 180 %(WoW 0)
    fileA = '2017-102n_O_A_BR_1to50_018perc_000deg_Test2.dat';
    fileB = '2017-102n_O_B_BR_1to50_018perc_000deg_Test4.dat';  
    fileD = '2017-102n_O_D_BR_1to50_018perc_000deg_Test6.dat'; 
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];     
end
dataA = fullfile(path,fileA);
dataB = fullfile(path,fileB);
dataD = fullfile(path,fileD);

% pressure data
psi2Pa  = 6894.76;
pdata_A = psi2Pa*importScanivalveDataFile(dataA);
pdata_B = psi2Pa*importScanivalveDataFile(dataB);
pdata_D = psi2Pa*importScanivalveDataFile(dataD);

% transfer function
path_TF = '/home/giacomol/Desktop/Research/windLoading/windTunnel/WoW/transferFunction';
file_TF = '2017-102n-Stanford-High-rise_FS_1652Hz_TS13_5in_plus_6in_CS12in_TFunc_Test11_1_10min_FF_250Hz_FSO_520Hz.mat';
TFFName = fullfile(path_TF,file_TF);

% apply transfer function
pdata_A = applyTF(pdata_A, TFFName, 520, 250);
pdata_B = applyTF(pdata_B, TFFName, 520, 250);
pdata_D = applyTF(pdata_D, TFFName, 520, 250);

% mean pressure from Z taps
pmeanZ_A = 0;%mean(nanmean(pdata_A(:,Ntaps_A+1:end)));
pmeanZ_B = 0;%mean(nanmean(pdata_B(:,Ntaps_B+1:end)));
pmeanZ_D = 0;%mean(nanmean(pdata_D(:,Ntaps_B+1:end)));

% taps names
[~,txt,~] = xlsread('/home/storage/WoW_data/TestDocuments/Drawings/2017-102n-Tubing-Table.xlsx', 'Tile A');
taps_A    = txt(2:end,2);
[~,txt,~] = xlsread('/home/storage/WoW_data/TestDocuments/Drawings/2017-102n-Tubing-Table.xlsx', 'Tile B');
taps_B    = txt(2:Ntaps_B+1,2); 
taps_D    = taps_B; 

for i = 1:Ntaps_B
    taps_D{i}(1) = 'D';
end

% remove corrupted data
% [~, ~, index] = setxor({'A1112'}, taps_A(1:Ntaps_A,1:5));
[~, ~, index] = setxor({'A1514'}, taps_A(1:Ntaps_A,:));

taps_A     = taps_A(index,:);
taps_A_ref = taps_A_ref(index,:);
coords_A   = coords_A(index,:);
pdata_A    = pdata_A(:,index);

[~, ~, index] = setxor({''}, taps_B(1:Ntaps_B,:));

taps_B     = taps_B(index,:);
taps_B_ref = taps_B_ref(index,:);
coords_B   = coords_B(index,:);
pdata_B    = pdata_B(:,index);

% convert to array
taps_A = cell2mat(taps_A);
taps_B = cell2mat(taps_B);
taps_D = cell2mat(taps_D);

% re-arrange data in reference order
[~, index_A_ref, index_A] = intersect(taps_A_ref, taps_A, 'rows', 'stable');
[~, index_B_ref, index_B] = intersect(taps_B_ref, taps_B, 'rows', 'stable'); 
[~, index_D_ref, index_D] = intersect(taps_D_ref, taps_D, 'rows', 'stable'); 

pdata_A = pdata_A(:,index_A_ref);
pdata_B = pdata_B(:,index_B_ref);
pdata_D = pdata_D(:,index_D_ref);

% reorder coordinates and convert to meters
coords_A = coords_A(index_A_ref,:)/1000;
coords_B = coords_B(index_B_ref,:)/1000;
coords_D = coords_D(index_D_ref,:)/1000;

% convert coordinates to CFD reference system
coords_A(:,1) = 1+coords_A(:,1);
coords_B(:,1) = 1+coords_B(:,1);

if angle == 180;
    coords_A(:,1) = 1-coords_A(:,1);
    coords_B(:,1) = 1-coords_B(:,1);
end
coords_D(:,2) = 2-coords_D(:,2);

%% stagnation pressure
fileS   = '2017-102n_O_D_BR_1to50_018perc_000deg_Test6.dat';
dataS   = fullfile(path, fileS);
pdata_S = psi2Pa*importScanivalveDataFile(dataS);

% qA = mean(pdata_S(:,index_D(203))); qB = qA; qD = qA;

%% Pressure mean and variance
% tile A
pmean_A = nanmean(pdata_A,1);
pstd_A  = nanstd(pdata_A,1);

% prssure coefficients
cp_A     = pdata_A/qA;
cpmean_A = (pmean_A-pmeanZ_A)/qA;
cpstd_A  = pstd_A/qA;

% tile B
pmean_B = nanmean(pdata_B,1);
pstd_B  = nanstd(pdata_B,1);

% prssure coefficients
cp_B     = pdata_B/qB;
cpmean_B = (pmean_B-pmeanZ_B)/qB;
cpstd_B  = pstd_B/qB;

% tile D
pmean_D = nanmean(pdata_D,1);
pstd_D  = nanstd(pdata_D,1);

% prssure coefficients
cp_D     = pdata_D/qD;
cpmean_D = (pmean_D-pmeanZ_D)/qD;
cpstd_D  = pstd_D/qD;

% define structures
% tileA
tileA.CI          = 'no';
tileA.taps        = taps_A(index_A_ref,:);
tileA.coords      = coords_A;
tileA.time        = linspace(0, size(pdata_A,1)/fsamp, size(pdata_A,1));
tileA.timeHistory = cp_A;
tileA.mean        = cpmean_A';
tileA.std         = cpstd_A';
tileA.U           = U;
tileA.q           = qA;
%tileA.areaAverage = areaAveragedPressure(tileA, panelA);
%tileA.peak        = pressurePeak(tileA, 6, 0.22, 'cook');
%tileB
tileB.CI          = 'no';
tileB.taps        = taps_B(index_B_ref,:);
tileB.coords      = coords_B;
tileB.time        = linspace(0, size(pdata_B,1)/fsamp, size(pdata_B,1));
tileB.timeHistory = cp_B;
tileB.mean        = cpmean_B';
tileB.std         = cpstd_B';
tileB.U           = U;
tileB.q           = qB;
%tileB.areaAverage = areaAveragedPressure(tileB, panelB);
%tileB.peak        = pressurePeak(tileB, 6, 0.22, 'cook');
%tileD
tileD.CI          = 'no';
tileD.taps        = taps_D(index_D,:);
tileD.coords      = coords_D;
tileD.time        = linspace(0, size(pdata_D,1)/fsamp, size(pdata_D,1));
tileD.timeHistory = cp_D;
tileD.mean        = cpmean_D';
tileD.std         = cpstd_D';
tileD.U           = U;
tileD.q           = qB;
% tileD.areaAverage = areaAveragedPressure(tileD, panelB);
%tileD.peak        = pressurePeak(tileD, 6, 0.22, 'cook');
