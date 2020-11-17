function [tileA, tileB, tileD] = pressurePoliMi(angle)

    % Function that read .h5 file from PoliMi experiment of tiles A-B-D
    % and returns:
        
        %structure containing:
        %      coords: coordinates of pressure taps               [Ntaps,2]
        %        taps: name and order of pressure taps            [Ntaps,1]
        % timeHistory: time series of pressure coefficients       [Nsamp,Ntaps]
        %        mean: mean of pressure coefficients              [Ntaps,1]
        %         std: root mean square of pressure coefficients  [Ntaps,1]        
        
%% taps list - reference
tapsList   = load('/home/giacomol/Desktop/Research/windLoading/windTunnel/matlabCode/ArupGUI3.0/Orifizi-v2.mat');
coords_B   = tapsList.prese{2};
coords_D   = tapsList.prese{4};
taps_A_ref = cell2mat(tapsList.nomi{1}');
taps_B_ref = cell2mat(tapsList.nomi{2}');
taps_D_ref = cell2mat(tapsList.nomi{4}');

%% data
fsamp   = 500;
Ntaps_A = 224;
Ntaps_B = 223;

pathAB  = '/home/storage/PoliMi_data/TilesA-B/DATA/hdf5';
pathCD  = '/home/storage/PoliMi_data/TilesC-D/DATA/hdf5';

if angle == 0
    fileA = '1028_SIW_OriArup_P7000_Ang+000.h5';
    fileB = '1107_SIW_OriArup_P7000_Ang+000.h5'; 
    fileD = '2056_SIW_TileD_FA_P7000_Ang+000.h5';
    % area-averaged panels
    panelA = [0.96 1.00 1.94 2.00];
    panelB = [0.96 1.00 0.97 1.03];
elseif angle == 10
    fileA = '1026_SIW_OriArup_P7000_Ang+010.h5';
    fileB = '1109_SIW_OriArup_P7000_Ang+010.h5';  
    fileD = '2050_SIW_TileD_FA_P7000_Ang+010.h5';
    % area-averaged panels
    panelA = [0.96 1.00 1.94 2.00];
    panelB = [0.96 1.00 0.97 1.03];    
elseif angle == 20
    fileA = '1024_SIW_OriArup_P7000_Ang+020.h5';
    fileB = '1111_SIW_OriArup_P7000_Ang+020.h5';    
    fileD = '2044_SIW_TileD_FA_P7000_Ang+020.h5';
    % area-averaged panels
    panelA = [0.96 1.00 1.94 2.00];
    panelB = [0.96 1.00 0.97 1.03];
elseif angle == 30
    fileA = '1022_SIW_OriArup_P7000_Ang+030.h5';
    fileB = '1113_SIW_OriArup_P7000_Ang+030.h5';    
    fileD = '2038_SIW_TileD_FA_P7000_Ang+030.h5';
    % area-averaged panels
    panelA = [0.96 1.00 1.94 2.00];
    panelB = [0.96 1.00 0.97 1.03];
elseif angle == 135
    fileA = '1021_SIW_OriArup_P7000_Ang+135.h5';
    fileB = '1114_SIW_OriArup_P7000_Ang+135.h5';
    fileD = '2516_SIW_SF_TileD_FA_P7000_Ang-150.h5';
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];      
elseif angle == 150
    fileA = '1020_SIW_OriArup_P7000_Ang+150.h5';
    fileB = '1115_SIW_OriArup_P7000_Ang+150.h5';
    fileD = '2516_SIW_SF_TileD_FA_P7000_Ang-150.h5';
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];        
elseif angle == 170
    fileA = '1018_SIW_OriArup_P7000_Ang+170.h5';
    fileB = '1117_SIW_OriArup_P7000_Ang+170.h5';
    fileD = '2014_SIW_TileD_FA_P7000_Ang+170.h5';
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];    
elseif angle == 180
    fileA = '1016_SIW_OriArup_P7000_Ang-180.h5';
    fileB = '1119_SIW_OriArup_P7000_Ang-180.h5';   
    fileD = '2008_SIW_TileD_FA_P7000_Ang-180.h5';    
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];
elseif angle == 190
    fileA = '1013_SIW_OriArup_P7000_Ang-170.h5';
    fileB = '1121_SIW_OriArup_P7000_Ang-170.h5';   
    fileD = '2501_SIW_SF_TileD_FA_P7000_Ang-170.h5';
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];
elseif angle == 200
    fileA = '1011_SIW_OriArup_P7000_Ang-160.h5';
    fileB = '1123_SIW_OriArup_P7000_Ang-160.h5';   
    fileD = '2510_SIW_SF_TileD_FA_P7000_Ang-160.h5';
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];   
elseif angle == 210
    fileA = '1009_SIW_OriArup_P7000_Ang-150.h5';
    fileB = '1125_SIW_OriArup_P7000_Ang-150.h5';   
    fileD = '2516_SIW_SF_TileD_FA_P7000_Ang-150.h5';
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];
elseif angle == 225
    fileA = '1006_SIW_OriArup_P7000_Ang-135.h5';
    fileB = '1126_SIW_OriArup_P7000_Ang-135.h5';   
    fileD = '2516_SIW_SF_TileD_FA_P7000_Ang-150.h5';
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];    
elseif angle == 350
    fileA = '1030_SIW_OriArup_P7000_Ang-010.h5';
    fileB = '1105_SIW_OriArup_P7000_Ang-010.h5';   
    fileD = '2062_SIW_TileD_FA_P7000_Ang-010.h5'; 
    % area-averaged panels
    panelA = [0.00 0.04 1.94 2.00];
    panelB = [0.00 0.04 0.97 1.03];  
end
dataA = fullfile(pathAB,fileA);
dataB = fullfile(pathAB,fileB);
dataD = fullfile(pathCD,fileD);

% check if file is correct
h5disp(dataA)
% h5disp(dataB)
% h5disp(dataD)

% dynamic pressure
q = mean(h5read(dataA,'/flowData/q')) / 0.8813^2; % convert to 2m height

% velocity
U = mean(h5read(dataA,'/flowData/U'));

% pressure data
pdata_A = h5read(dataA,'/pressureData/pressureData');
pdata_B = h5read(dataB,'/pressureData/pressureData');
pdata_D = h5read(dataD,'/pressureData/pressureData');

% mean pressure from Z taps
pmeanZ_A = 0;%mean(nanmean(pdata_A(:,Ntaps_A+1:end)));
pmeanZ_B = 0;%mean(nanmean(pdata_B(:,Ntaps_B+1:end)));
pmeanZ_D = 0;%mean(nanmean(pdata_D(:,Ntaps_B+1:end)));

% taps names
taps_A = cell2mat(h5read(dataA,'/metadata/tapsNames'));
taps_B = cell2mat(h5read(dataB,'/metadata/tapsNames')); 
taps_D = cell2mat(h5read(dataD,'/metadata/tapsNames')); 

% taps coordinates
coords_A = h5read(dataA,'/metadata/tapsCoords');
% coords_D = h5read(dataD,'/metadata/tapsCoords');

% remove corrupted data
% [~, ~, index] = setxor({'A1112'}, taps_A(1:Ntaps_A,1:5));
[~, ~, index] = setxor({'A1514'}, taps_A(1:Ntaps_A,1:5));

taps_A     = taps_A(index,:);
taps_A_ref = taps_A_ref(index,:);
coords_A   = coords_A(index,:);
pdata_A    = pdata_A(:,index);

% re-arrange data in reference order
[~, index_A_ref, index_A] = intersect(taps_A_ref, taps_A(:,1:5), 'rows', 'stable');
[~, index_B_ref, index_B] = intersect(taps_B_ref, taps_B(:,1:5), 'rows', 'stable'); 
[~, index_D_ref, index_D] = intersect(taps_D_ref, taps_D(:,1:5), 'rows', 'stable'); 

pdata_A = pdata_A(:,index_A_ref);
pdata_B = pdata_B(:,index_B);
pdata_D = pdata_D(:,index_D);

% reorder coordinates and convert to meters
coords_A = coords_A(index_A_ref,:)/1000;
coords_B = coords_B/1000;
coords_D = coords_D/1000;

% convert coordinates to CFD reference system
coords_A(:,1) = 1+coords_A(:,1);
coords_B(:,1) = 1+coords_B(:,1);
coords_D(:,2) = 2-coords_D(:,2);

if angle > 90 && angle < 270;
    coords_A(:,1) = 1-coords_A(:,1);
    coords_B(:,1) = 1-coords_B(:,1);
    coords_D(:,1) = 0.3-coords_D(:,1);
end

%% stagnation pressure
fileS   = '2008_SIW_TileD_FA_P7000_Ang-180.h5';
dataS   = fullfile(pathCD, fileS);
pdata_S = h5read(dataS,'/pressureData/pressureData');

% qA = mean(pdata_S(:,index_D(203))); qB = qA; qD = qA;

%% Pressure mean and root mean square
% tile A
pmean_A = nanmean(pdata_A,1);
pstd_A  = nanstd(pdata_A,1);

% prssure coefficients
cp_A     = pdata_A/q;
cpmean_A = (pmean_A-pmeanZ_A)/q;
cpstd_A  = pstd_A/q;

% tile B
pmean_B = nanmean(pdata_B,1);
pstd_B  = nanstd(pdata_B,1);

% prssure coefficients
cp_B     = pdata_B/q;
cpmean_B = (pmean_B-pmeanZ_B)/q;
cpstd_B  = pstd_B/q;

% tile D
pmean_D = nanmean(pdata_D,1);
pstd_D  = nanstd(pdata_D,1);

% prssure coefficients
cp_D     = pdata_D/q;
cpmean_D = (pmean_D-pmeanZ_D)/q;
cpstd_D  = pstd_D/q;

% define structures
% tileA
tileA.CI          = 'no';
tileA.taps        = taps_A(index_A_ref,1:5);
tileA.coords      = coords_A;
tileA.time        = linspace(0, size(pdata_A,1)/fsamp, size(pdata_A,1));
tileA.timeHistory = cp_A;
tileA.mean        = cpmean_A';
tileA.std         = cpstd_A';
tileA.U           = U;
tileA.q           = q;
%tileA.areaAverage = areaAveragedPressure(tileA, panelA);
%tileA.peak        = pressurePeak(tileA, 1, 0.22, 'cook');

% [tileA.labels, tileA.panels] = designPressureLabels('on', tileA);

%tileB
tileB.CI          = 'no';
tileB.taps        = taps_B(index_B,1:5);
tileB.coords      = coords_B;
tileB.time        = linspace(0, size(pdata_B,1)/fsamp, size(pdata_B,1));
tileB.timeHistory = cp_B;
tileB.mean        = cpmean_B';
tileB.std         = cpstd_B';
tileB.U           = U;
tileB.q           = q;
%tileB.areaAverage = areaAveragedPressure(tileB, panelB);
%tileB.peak        = pressurePeak(tileB, 1, 0.22, 'cook');

% [tileB.labels, tileB.panels] = designPressureLabels('on', tileB);

%tileD
tileD.CI          = 'no';
tileD.taps        = taps_D(index_D,1:5);
tileD.coords      = coords_D;
tileD.time        = linspace(0, size(pdata_D,1)/fsamp, size(pdata_D,1));
tileD.timeHistory = cp_D;
tileD.mean        = cpmean_D';
tileD.std         = cpstd_D';
tileD.U           = U;
tileD.q           = q;
% tileD.areaAverage = areaAveragedPressure(tileD);
%tileD.peak        = pressurePeak(tileD, 1, 0.22, 'cook');

% [tileD.labels, tileD.panels] = designPressureLabels('on', tileD);
