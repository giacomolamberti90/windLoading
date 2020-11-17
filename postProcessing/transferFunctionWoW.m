%% Apply transfer function to WoW data

clear; close all; clc

angle = 0;

%% raw data
path_raw = '/home/storage/WoW_data/Raw_Data_Untouched/Scanivalve_Data/Braced';

if angle == 0 %(WoW 180)
    fileA = '2017-102n_O_A_BR_1to50_018perc_180deg_Test2.dat';
    fileB = '2017-102n_O_B_BR_1to50_018perc_180deg_Test4.dat'; 
    fileD = '2017-102n_O_D_BR_1to50_018perc_180deg_Test6.dat';
elseif angle == 10 %(WoW 170)
    fileA = '2017-102n_O_A_BR_1to50_018perc_170deg_Test2.dat';
    fileB = '2017-102n_O_B_BR_1to50_018perc_170deg_Test4.dat';  
    fileD = '2017-102n_O_D_BR_1to50_018perc_170deg_Test6.dat';
elseif angle == 20 %(WoW 160)
    fileA = '2017-102n_O_A_BR_1to50_018perc_160deg_Test2.dat';
    fileB = '2017-102n_O_B_BR_1to50_018perc_160deg_Test4.dat';    
    fileD = '2017-102n_O_D_BR_1to50_018perc_160deg_Test6.dat';
elseif angle == 180 %(WoW 0)
    fileA = '2017-102n_O_A_BR_1to50_018perc_000deg_Test2.dat';
    fileB = '2017-102n_O_B_BR_1to50_018perc_000deg_Test4.dat';  
    fileD = '2017-102n_O_D_BR_1to50_018perc_000deg_Test6.dat';
end
dataA = fullfile(path_raw, fileA);
dataB = fullfile(path_raw, fileB);
dataD = fullfile(path_raw, fileD);

% pressure data
psi2Pa  = 6894.76;
pdata_A = importScanivalveDataFile(dataA);
pdata_B = importScanivalveDataFile(dataB);
pdata_D = importScanivalveDataFile(dataD);

%% transfer function
path_TF   = '/home/giacomol/Desktop/Research/windLoading/windTunnel/WoW/transferFunction';
file_TF_B = '2017-102n-Stanford-High-rise_FS_1652Hz_TS13_5in_plus_6in_CS12in_TFunc_Test11_15min_FF_250Hz_FSO_520Hz.mat';
file_TF_A = '2017-102n-Stanford-High-rise_FS_1652Hz_TS77_25in_plus_6in_CS12in_TFunc_Test12_15min_FF_250Hz_FSO_520Hz.mat';
TFFName_A = fullfile(path_TF,file_TF_A);
TFFName_B = fullfile(path_TF,file_TF_B);

%% apply transfer function
pdata_TF_A = applyTF(pdata_A, TFFName_A);
pdata_TF_B = applyTF(pdata_B, TFFName_B);
pdata_TF_D = applyTF(pdata_D, TFFName_B);

pdata_A = psi2Pa*pdata_A;
pdata_B = psi2Pa*pdata_B;
pdata_D = psi2Pa*pdata_D;

pdata_TF_A = psi2Pa*pdata_TF_A;
pdata_TF_B = psi2Pa*pdata_TF_B;
pdata_TF_D = psi2Pa*pdata_TF_D;

figure
hold on
plot(mean(pdata_A),'.')
plot(mean(pdata_TF_A),'.')

figure
hold on
plot(std(pdata_A),'.')
plot(std(pdata_TF_A),'.')

figure
hold on
plot(mean(pdata_B),'.')
plot(mean(pdata_TF_B),'.')

figure
hold on
plot(std(pdata_B),'.')
plot(std(pdata_TF_B),'.')

%% write filtered data
% path_post = '/home/storage/WoW_data/Post_Processed_Data/Scanivalve_Data/Braced';
% 
% dataA = fullfile(path_post, fileA);
% dataB = fullfile(path_post, fileB);
% dataD = fullfile(path_post, fileD);
% 
% % pressure data
% pdata_TF_A = 6894.76*load(dataA);
% pdata_TF_B = 6894.76*load(dataB);
% pdata_TF_D = 6894.76*load(dataD);

