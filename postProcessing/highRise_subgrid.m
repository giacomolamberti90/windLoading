%% POST-PROCESSING PRESSURE ON HIGH-RISE ----------------------------------

clear; close all; clc

%% Load data --------------------------------------------------------------

%% Wind Tunnel ------------------------------------------------------------
% PoliMi
% 0-180 (left/right face in CFD)
[tileA_poli_0,   tileB_poli_0,   tileD_poli_0]   = pressurePoliMi(0);
[tileA_poli_20,  tileB_poli_20,  tileD_poli_20]  = pressurePoliMi(20);
[tileA_poli_180, tileB_poli_180, tileD_poli_180] = pressurePoliMi(180);

%% LES convergence --------------------------------------------------------
% path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/20deg/coarse/dynamicK/workdir.14';
% path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/20deg/coarse/dynamicLagrangian/workdir.14';
% path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/20deg/coarse/smagorinsky/workdir.14';

% convergenceLES(path_LES, 20, 1, 1)
% convergenceLES(path_LES, 180, 1, 1)

%% LES --------------------------------------------------------------------
% Filter LES data using experiment time resolution
deltaT = 0.0001;
deltaT_exp = deltaT;

% 00deg
% dynamicK
path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/00deg/coarse/dynamicK/workdir.14';
[tileA_LESk_0, tileB_LESk_0, ~] = pressureLES(path_LES, 0, deltaT, deltaT_exp);
[tileA_LESk_180, tileB_LESk_180, ~] = pressureLES(path_LES, 180, deltaT, deltaT_exp);

% dynamiclagrangian
path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/00deg/coarse/dynamicLagrangian/workdir.14';
[tileA_LESl_0, tileB_LESl_0, ~] = pressureLES(path_LES, 0, deltaT, deltaT_exp);
[tileA_LESl_180, tileB_LESl_180, ~] = pressureLES(path_LES, 180, deltaT, deltaT_exp);

% smagorinsky
path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/00deg/coarse/smagorinsky/workdir.14';
[tileA_LESs_0, tileB_LESs_0, ~] = pressureLES(path_LES, 0, deltaT, deltaT_exp);
[tileA_LESs_180, tileB_LESs_180, ~] = pressureLES(path_LES, 180, deltaT, deltaT_exp);

tileA_LESk = combineTilesA(tileA_LESk_0, tileA_LESk_180);
tileA_LESl = combineTilesA(tileA_LESl_0, tileA_LESl_180);
tileA_LESs = combineTilesA(tileA_LESs_0, tileA_LESs_180);

% 20deg
% dynamicK
path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/20deg/coarse/dynamicK/workdir.14';
[tileA_LESk_20, tileB_LESk_20, ~] = pressureLES(path_LES, 20, deltaT, deltaT_exp);

% dynamiclagrangian
path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/20deg/coarse/dynamicLagrangian/workdir.14';
[tileA_LESl_20, tileB_LESl_20, ~] = pressureLES(path_LES, 20, deltaT, deltaT_exp);

% smagorinsky
path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/20deg/coarse/smagorinsky/workdir.14';
[tileA_LESs_20, tileB_LESs_20, ~] = pressureLES(path_LES, 20, deltaT, deltaT_exp);

% Charles
path_LES = '/home/giacomol/Desktop/Research/windLoading/windTunnel/charles/';
[tileA_LESc_0, tileB_LESc_0] = pressureCharlesLES(path_LES, 0);
[tileA_LESc_20, tileB_LESc_20] = pressureCharlesLES(path_LES, 20);
[tileA_LESc_180, tileB_LESc_180] = pressureCharlesLES(path_LES, 180);

tileA_LESc = combineTilesA(tileA_LESc_0, tileA_LESc_180);

return

%% Pressure contours ------------------------------------------------------
% select plot
stat = 'peak'; % select variable to plot
ax = [-3 0];
path = '/home/giacomol/Desktop/Research/slides/02-06-20/figures';

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESc_20, stat, ax);
hold on
pressureContour('no', tileB_LESc_20, stat, ax);
% pressureContour('no', tileA_LESc_180, stat, ax);
% pressureContour('no', tileB_LESc_180, stat, ax);
title('$$Charles$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(path, 'cppeak_charles_20deg_contour.png'))

%%
% 00deg -------------------------------------------------------------------
figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_poli_0, stat, ax);
hold on
pressureContour('no', tileB_poli_0, stat, ax);
pressureContour('no', tileA_poli_180, stat, ax);
pressureContour('no', tileB_poli_180, stat, ax);
title('$$PoliMi$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(path, 'cppeak_poli_00deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESk_0, stat, ax);
hold on
pressureContour('no', tileB_LESk_0, stat, ax);
pressureContour('no', tileA_LESk_180, stat, ax);
pressureContour('no', tileB_LESk_180, stat, ax);
title('$$Dynamic K$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(path, 'cppeak_dynK_00deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESl_0, stat, ax);
hold on
pressureContour('no', tileB_LESl_0, stat, ax);
pressureContour('no', tileA_LESl_180, stat, ax);
pressureContour('no', tileB_LESl_180, stat, ax);
title('$$Dynamic Lagrangian$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(path, 'cppeak_dynL_00deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESs_0, stat, ax);
hold on
pressureContour('no', tileB_LESs_0, stat, ax);
pressureContour('no', tileA_LESs_180, stat, ax);
pressureContour('no', tileB_LESs_180, stat, ax);
title('$$Smagorinsky$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(path, 'cppeak_smag_00deg_contour.png'))

%%
stat = 'peak'; % select variable to plot
ax = [-5 0];

% 20deg -------------------------------------------------------------------
figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_poli_20, stat, ax);
pressureContour('no', tileB_poli_20, stat, ax);
title('$$PoliMi$$','interpreter','latex', 'fontsize', 22)
% axis([0.8 1 1.8 2])
% saveas(gcf, fullfile(path, 'cprms_poli_20deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESk_20, stat, ax);
% pressureContour('no', tileB_LESk_20, stat, ax);
title('$$Dynamic K$$','interpreter','latex', 'fontsize', 22)
% axis([0.8 1 1.8 2])
% saveas(gcf, fullfile(path, 'cprms_dynK_20deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESl_20, stat, ax);
% pressureContour('no', tileB_LESl_20, stat, ax);
title('$$Dynamic Lagrangian$$','interpreter','latex', 'fontsize', 22)
% axis([0.8 1 1.8 2])
% saveas(gcf, fullfile(path, 'cprms_dynL_20deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESs_20, stat, ax);
% pressureContour('no', tileB_LESs_20, stat, ax);
title('$$Smagorinsky$$','interpreter','latex', 'fontsize', 22)
% axis([0.8 1 1.8 2])
% saveas(gcf, fullfile(path, 'cprms_smag_20deg_contour.png'))

%% Pressure profiles ------------------------------------------------------
stat = 'peak';

% Tile A
figure('rend','painters','pos', [0 0 550 420]);
hold on

pressureProfile_all(tileA_LESs, 'A', stat, '-ok');
pressureProfile_all(tileA_LESk, 'A', stat, '-ob');
pressureProfile_all(tileA_LESl, 'A', stat, '-oc');
pressureProfile_all(tileA_LESc, 'A', stat, '--sm');

% pressureProfile(tileA_LESs_20, stat, '-ok'); 
% pressureProfile(tileA_LESk_20, stat, '-ob'); 
% pressureProfile(tileA_LESl_20, stat, '-oc'); 
% pressureProfile(tileA_LESc_20, stat, '--sm');

pressureProfile(tileA_poli_0, stat, '.r');
pressureProfile(tileA_poli_180, stat, '.r');
ylim([-3 0])

% legend('Smagorinsky', 'DynamicK', 'Lagrangian', 'Charles', 'PoliMi'); legend boxoff
% saveas(gcf, fullfile(path, 'cppeak_A_00deg_profiles.png'))

% Tile B
figure('rend','painters','pos', [0 0 550 420]);
hold on

pressureProfile(tileB_LESs_0, stat, '-ok'); 
pressureProfile(tileB_LESk_0, stat, '-ob'); 
pressureProfile(tileB_LESl_0, stat, '-oc'); 
pressureProfile(tileB_LESc_0, stat, '--sm');

pressureProfile(tileB_poli_0, stat, '.r');

pressureProfile(tileB_LESs_180, stat, '-ok'); 
pressureProfile(tileB_LESk_180, stat, '-ob'); 
pressureProfile(tileB_LESl_180, stat, '-oc'); 
pressureProfile(tileB_LESc_180, stat, '--sm'); 

pressureProfile(tileB_poli_180, stat, '.r');

ylim([-3 0])
legend('Smagorinsky', 'DynamicK', 'Lagrangian', 'Charles', 'PoliMi'); legend boxoff
% saveas(gcf, fullfile(path, 'cppeak_B_00deg_profiles.png'))

%% Pressure time-series ---------------------------------------------------
list = 'A0301';

tin = 50;
tend = 52;

% 00deg -------------------------------------------------------------------
figure('rend','painters','pos', [0 0 550 420]);
pressureTimeHistory('no', tileA_poli_0, 'r', -5, [tin tend]);
pressureTimeHistory('no', tileA_LESk_0, 'b', -5, [tin tend]);
pressureTimeHistory('no', tileA_LESl_0, 'c', -5, [tin tend]);
pressureTimeHistory('no', tileA_LESs_0, 'k', -5, [tin tend]);
% saveas(gcf, fullfile(path, 'cp_A0301_180deg.png'))

list = 'B0810';

figure('rend','painters','pos', [0 0 550 420]);
pressureTimeHistory('no', tileB_poli_0, 'r', -5, [tin tend]);
pressureTimeHistory('no', tileB_LESk_0, 'b', -5, [tin tend]);
pressureTimeHistory('no', tileB_LESl_0, 'c', -5, [tin tend]);
pressureTimeHistory('no', tileB_LESs_0, 'k', -5, [tin tend]);
% legend('PoliMi','dynamick', 'dynamicLagrangian', 'Smagorinsky') 
% legend boxoff
% saveas(gcf, fullfile(path, 'cp_B0810_180deg.png'))


%% Pressure spectra -------------------------------------------------------
% 00deg
list = 'A0205';

LES = [5,14,23];

figure
pressurePSD(tileA_poli_0, 500, 'adi', 'r');
hold on
pressurePSD(tileA_LESs_0, 1/0.001, 'adi', 'k');
pressurePSD(tileA_LESk_0, 1/0.001, 'adi', 'b');
pressurePSD(tileA_LESl_0, 1/0.0002, 'adi', 'c');

legend('PoliMi', 'Smagorinsky', 'DynamicK', 'Lagrangian'); legend boxoff

list = 'B0810';

figure
pressurePSD(tileB_poli_0, 500, 'adi', 'r');
hold on
pressurePSD(tileB_LESs_0, 1/0.001, 'adi', 'k');
pressurePSD(tileB_LESk_0, 1/0.001, 'adi', 'b');
pressurePSD(tileB_LESl_0, 1/0.0002, 'adi', 'c');

%%
% 20deg
list = 'A0205';

tin = 4;
tend = 6;

X_values = unique(tileA_LESk_20.coords(:,1));
Y_values = unique(tileA_LESk_20.coords(:,2));
X_values = [0.0156; 0.0944; 0.5088; 0.9965];

for i = 3%:length(X_values)

    [row, col] = find((tileA_LESk_20.coords(:,2) == Y_values(end-20)) & ...
                       round(tileA_LESk_20.coords(:,1), 3) == round(X_values(i), 3));

    %figure('rend','painters','pos', [0 0 550 420]);
    % pressureTimeHistory('no', tileA_poli_20, 'r', -15, [tin tend], list);
    pressureTimeHistory('no', tileA_LESk_20, 'b', -7, [tin tend], list, tileA_LESk_20.timeHistory(:,row));
    pressureTimeHistory('no', tileA_LESl_20, 'c', -7, [tin tend], list, tileA_LESl_20.timeHistory(:,row));
    pressureTimeHistory('no', tileA_LESs_20, 'k', -7, [tin tend], list, tileA_LESs_20.timeHistory(:,row));
    ylim([-1.5 0])
    % legend('DynamicK', 'Lagrangian', Smagorinsky); legend boxoff
%     saveas(gcf, fullfile(path, sprintf('cp_X%i_20deg.png', i)))

    % figure
    % % pressurePSD(tileA_poli_20, 500, 'adi', 'r');
    % pressurePSD(tileA_LESs_20, 1/0.0005, 'adi', 'k');
    % hold on
    % pressurePSD(tileA_LESk_20, 1/0.001, 'adi', 'b');
    % pressurePSD(tileA_LESl_20, 1/0.0004, 'adi', 'c');
    % % legend('Smagorinsky', 'DynamicK', 'Lagrangian'); legend boxoff

end

%% EVA vs area ------------------------------------------------------------
POE = 0.22;
ang = 0;

% PoliMi
[tileA_poli_0.design, area] = designPressure('no', tileA_poli_0, ang, 6, POE);

% LES
[tileA_LESk_0.design, ~] = designPressure('no', tileA_LESk_0, ang, 6, POE);
[tileA_LESl_0.design, ~] = designPressure('no', tileA_LESl_0, ang, 6, POE);
[tileA_LESs_0.design, ~] = designPressure('no', tileA_LESs_0, ang, 6, POE);

figure
semilogx(area, tileA_poli_0.design, '.r', 'markersize', 15)
hold on
semilogx(area, tileA_LESs_0.design, '-ok')
semilogx(area, tileA_LESk_0.design, '-ob')
semilogx(area, tileA_LESl_0.design, '-oc')
legend('PoliMi', 'Smagorinsky', 'DynamicK', 'Lagrangian'); legend boxoff
set(gca,'fontsize',22)
ylim([-1 0])

%%
ang = 20;

% PoliMi
[tileA_poli_20.design, area] = designPressure('no', tileA_poli_20, ang, 6, POE);

% LES
[tileA_LESk_20.design, ~] = designPressure('no', tileA_LESk_20, ang, 6, POE);
[tileA_LESl_20.design, ~] = designPressure('no', tileA_LESl_20, ang, 6, POE);
[tileA_LESs_20.design, ~] = designPressure('no', tileA_LESs_20, ang, 6, POE);

figure
semilogx(area, tileA_poli_20.design, '.r', 'markersize', 15)
hold on
semilogx(area, tileA_LESs_20.design, '-ok')
semilogx(area, tileA_LESk_20.design, '-ob')
semilogx(area, tileA_LESl_20.design, '-oc')
legend('PoliMi', 'Smagorinsky', 'DynamicK', 'Lagrangian'); legend boxoff
set(gca,'fontsize',22)
ylim([-3 0])
