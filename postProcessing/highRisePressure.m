%% POST-PROCESSING PRESSURE ON HIGH-RISE ----------------------------------

clear; close all; clc

%% Load data ------------------------------------- -------------------------

%% Wind Tunnel ------------------------------------------------------------
% PoliMi
% 0-180 (left/right face in CFD)
[tileA_poli_0,   tileB_poli_0,   tileD_poli_0]   = pressurePoliMi(0);
[tileA_poli_20,  tileB_poli_20,  tileD_poli_20]  = pressurePoliMi(20);
[tileA_poli_180, tileB_poli_180, tileD_poli_180] = pressurePoliMi(180);

return

% 10-170 (left face in CFD)
[tileA_poli_10,  tileB_poli_10,  tileD_poli_10]  = pressurePoliMi(10);
[tileA_poli_170, tileB_poli_170, tileD_poli_170] = pressurePoliMi(170);

% 20 (left face in CFD) 160 N/A
[tileA_poli_20, tileB_poli_20, tileD_poli_20] = pressurePoliMi(20);

% 30-150 (left face in CFD)
[tileA_poli_30,  tileB_poli_30,  tileD_poli_30]  = pressurePoliMi(30);
[tileA_poli_150, tileB_poli_150, tileD_poli_150] = pressurePoliMi(150);

% 190-350 (right face in CFD at 10deg)
[tileA_poli_190, tileB_poli_190, tileD_poli_190] = pressurePoliMi(190);
[tileA_poli_350, tileB_poli_350, tileD_poli_350] = pressurePoliMi(350);

% 190-350 (right face in CFD at 20deg) 340 N/A
[tileA_poli_200, tileB_poli_200, tileD_poli_200] = pressurePoliMi(200);

% 210 (right face in CFD at 30deg) 330 N/A
[tileA_poli_200, tileB_poli_200, tileD_poli_200] = pressurePoliMi(200);

%% RANS UQ ----------------------------------------------------------------
highRise_base_0 = pressureUQ_all(0, 'bas');
highRise_c1c_0  = pressureUQ_all(0, 'c1c');
highRise_c2c_0  = pressureUQ_all(0, 'c2c');
highRise_c3c_0  = pressureUQ_all(0, 'c3c');

highRise_base_20 = pressureUQ_all(20, 'bas');
highRise_c1c_20  = pressureUQ_all(20, 'c1c');
highRise_c2c_20  = pressureUQ_all(20, 'c2c');
highRise_c3c_20  = pressureUQ_all(20, 'c3c');

highRise_comb_0  = pressureCombinedUQ_all(highRise_base_0,  highRise_c1c_delta0p1_0,  highRise_c2c_delta0p5_0,  highRise_c3c_delta0p2_0);
highRise_comb_20 = pressureCombinedUQ_all(highRise_base_20, highRise_c1c_20, highRise_c2c_20, highRise_c3c_20);

return

% % RANS - aleatoric UQ
% [tileA_UQ_0, tileB_UQ_0]     = pressureUQ(0, 'mean');
% [tileA_UQ_20,  tileB_UQ_20]  = pressureUQ(20, 'mean');
% [tileA_UQ_180, tileB_UQ_180] = pressureUQ(180, 'mean');
% 
% % RANS - epistemic UQ - c1c
% [tileA_c1c_0,   tileB_c1c_0]   = pressureEpistemicUQ(0, 'c1c', 'mean');
% [tileA_c1c_20,  tileB_c1c_20]  = pressureEpistemicUQ(20, 'c1c', 'mean');
% [tileA_c1c_180, tileB_c1c_180] = pressureEpistemicUQ(180, 'c1c', 'mean');
% 
% % RANS - epistemic UQ - c2c
% [tileA_c2c_0,   tileB_c2c_0]   = pressureEpistemicUQ(0, 'c2c', 'mean');
% [tileA_c2c_05_20,  tileB_c2c_05_20]  = pressureEpistemicUQ(20, 'c2c', 'mean');
% [tileA_c2c_180, tileB_c2c_180] = pressureEpistemicUQ(180, 'c2c', 'mean');
% 
% % RANS - epistemic UQ - c3c
% [tileA_c3c_0,   tileB_c3c_0]   = pressureEpistemicUQ(0, 'c3c', 'mean');
% [tileA_c3c_20,  tileB_c3c_20]  = pressureEpistemicUQ(20, 'c3c', 'mean');
% [tileA_c3c_180, tileB_c3c_180] = pressureEpistemicUQ(180, 'c3c', 'mean');
% 
% [tileA_comb_0,   tileB_comb_0]  = pressureCombinedUQ(tileA_UQ_0, tileA_c1c_0, tileA_c2c_0, tileA_c3c_0,...
%                                                      tileB_UQ_0, tileB_c1c_0, tileB_c2c_0, tileB_c3c_0);
% 
% [tileA_comb_20,  tileB_comb_20] = pressureCombinedUQ(tileA_UQ_20, tileA_c1c_20, tileA_c2c_20, tileA_c3c_20,...
%                                                      tileB_UQ_20, tileB_c1c_20, tileB_c2c_20, tileB_c3c_20);
%                                                  
% [tileA_comb_180, tileB_comb_180] = pressureCombinedUQ(tileA_UQ_180, tileA_c1c_180, tileA_c2c_180, tileA_c3c_180,...
%                                                       tileB_UQ_180, tileB_c1c_180, tileB_c2c_180, tileB_c3c_180);

%% LES --------------------------------------------------------------------
path_SENS = '/home/storage/LES/highRise/fine/00deg/Smagorinsky/';
% path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/80deg';
% convergenceLES(path_LES, 80, 1, 1)
% 
% return
% 
% for i = 1:27
%     sprintf('workdir.%i', i)
%     path_LES  = fullfile(path_SENS, sprintf('workdir.%i', i));
%     convergenceLES(path_LES, 20, 1, 1)
%     %convergenceLES(path_LES, 180, 1, 1)
% end
% 
% return
% Filter LES data using experiment time resolution
deltaT = 0.0001;
deltaT_exp = deltaT; %tileA_poli_0.time(2);

for i = 14%1:27
    
    sprintf('workdir = %i', i)
    path_LES = fullfile(path_SENS, sprintf('workdir.%i', i));

    tic
    [tileA_LES_0{i}, tileB_LES_0{i}, ~] = pressureLES(path_LES, 0, deltaT, deltaT_exp);
    [tileA_LES_180{i}, tileB_LES_180{i}, ~] = pressureLES(path_LES, 180, deltaT, deltaT_exp);
    toc
end

return

path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/00deg/coarse/dynamicK/workdir.14';
[tileA_coarse_0, tileB_coarse_0, ~] = pressureLES(path_LES, 0, deltaT, deltaT_exp);
[tileA_coarse_180, tileB_coarse_180, ~] = pressureLES(path_LES, 180, deltaT, deltaT_exp);

path_LES = fullfile(path_SENS, sprintf('workdir.%i', 14));
path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/00deg/fine/dynamicK/workdir.14';
[tileA_fine_0, tileB_fine_0, ~] = pressureLES(path_LES, 0, deltaT, deltaT_exp);
[tileA_fine_180, tileB_fine_180, ~] = pressureLES(path_LES, 180, deltaT, deltaT_exp);

path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/20deg/coarse/dynamicK/workdir.14';
[tileA_LESk_20, tileB_LESk_20, ~] = pressureLES(path_LES, 20, deltaT, deltaT_exp);

path_LES = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/45deg/coarse/dynamicK/workdir.14';
[tileA_LESk_45, tileB_LESk_45, ~] = pressureLES(path_LES, 45, deltaT, deltaT_exp);

%% pressure peaks
for i = 1:27
    
    sprintf('workdir = %i', i)    
    %tileA_LES_0{i}.peak = pressurePeak(tileA_LES_0{i}, window, .22, method);
    %tileB_LES_0{i}.peak = pressurePeak(tileB_LES_0{i}, window, .22, method);
    
    %tileA_LES_180{i}.peak = pressurePeak(tileA_LES_180{i}, window, .22, method);
    %tileB_LES_180{i}.peak = pressurePeak(tileB_LES_180{i}, window, .22, method);

    tileA{i} = combineTilesA(tileA_LES_0{i}, tileA_LES_180{i});
    tileB{i} = combineTilesA(tileB_LES_0{i}, tileB_LES_180{i});
end

[tileA_sens_0,   tileB_sens_0,   tileD_sens_0]   = pressureCombinedLES(tileA_LES_0, tileB_LES_0, tileB_LES_0);
[tileA_sens_180, tileB_sens_180, tileD_sens_180] = pressureCombinedLES(tileA_LES_180, tileB_LES_180, tileB_LES_180);

tileA_sens = combineTilesA(tileA_sens_0, tileA_sens_180);
tileB_sens = combineTilesA(tileB_sens_0, tileB_sens_180);

%% compute peaks
window = 6;
POE = 0.22;
method = 'cook';

tileA_1_0 = tileA_poli_0; tileA_1_180 = tileA_poli_180;
tileB_1_0 = tileB_poli_0; tileB_1_180 = tileB_poli_180;

window = 1;
tileA_1_0.peak = pressurePeak(tileA_1_0, window, POE, method)
tileB_1_0.peak = pressurePeak(tileB_1_0, window, POE, method)
tileA_1_180.peak = pressurePeak(tileA_1_180, window, POE, method)
tileB_1_180.peak = pressurePeak(tileB_1_180, window, POE, method)

tileA_2_0 = tileA_poli_0; tileA_2_180 = tileA_poli_180;
tileB_2_0 = tileB_poli_0; tileB_2_180 = tileB_poli_180;

window = 2;
tileA_2_0.peak = pressurePeak(tileA_2_0, window, POE, method)
tileB_2_0.peak = pressurePeak(tileB_2_0, window, POE, method)
tileA_2_180.peak = pressurePeak(tileA_2_180, window, POE, method)
tileB_2_180.peak = pressurePeak(tileB_2_180, window, POE, method)

tileA_3_0 = tileA_poli_0; tileA_3_180 = tileA_poli_180;
tileB_3_0 = tileB_poli_0; tileB_3_180 = tileB_poli_180;

window = 3;
tileA_3_0.peak = pressurePeak(tileA_3_0, window, POE, method)
tileB_3_0.peak = pressurePeak(tileB_3_0, window, POE, method)
tileA_3_180.peak = pressurePeak(tileA_3_180, window, POE, method)
tileB_3_180.peak = pressurePeak(tileB_3_180, window, POE, method)

tileA_4_0 = tileA_poli_0; tileA_4_180 = tileA_poli_180;
tileB_4_0 = tileB_poli_0; tileB_4_180 = tileB_poli_180;

window = 4;
tileA_4_0.peak = pressurePeak(tileA_4_0, window, POE, method)
tileB_4_0.peak = pressurePeak(tileB_4_0, window, POE, method)
tileA_4_180.peak = pressurePeak(tileA_4_180, window, POE, method)
tileB_4_180.peak = pressurePeak(tileB_4_180, window, POE, method)

tileA_5_0 = tileA_poli_0; tileA_5_180 = tileA_poli_180;
tileB_5_0 = tileB_poli_0; tileB_5_180 = tileB_poli_180;

window = 5;
tileA_5_0.peak = pressurePeak(tileA_5_0, window, POE, method)
tileB_5_0.peak = pressurePeak(tileB_5_0, window, POE, method)
tileA_5_180.peak = pressurePeak(tileA_5_180, window, POE, method)
tileB_5_180.peak = pressurePeak(tileB_5_180, window, POE, method)

stat = 'peak';
ax = [-3, 0];

% 00deg -------------------------------------------------------------------
% figure('rend','painters','pos', [0 0 400 400]);
% pressureContour('no', tileA_1_0, stat, ax);
% hold on
% pressureContour('no', tileB_1_0, stat, ax);
% pressureContour('no', tileA_1_180, stat, ax);
% pressureContour('no', tileB_1_180, stat, ax);
% title('$$Coarse$$','interpreter','latex', 'fontsize', 22)
% % saveas(gcf, fullfile(path, sprintf('cp_%s_%s_fine_00deg_contour.png', stat, s)))
% 
% figure('rend','painters','pos', [0 0 400 400]);
% pressureContour('no', tileA_2_0, stat, ax);
% hold on
% pressureContour('no', tileB_2_0, stat, ax);
% pressureContour('no', tileA_2_180, stat, ax);
% pressureContour('no', tileB_2_180, stat, ax);
% title('$$Fine$$','interpreter','latex', 'fontsize', 22)
% % saveas(gcf, fullfile(path, sprintf('cp_%s_%s_coarse_00deg_contour.png', stat, s)))

% Tile A
figure('rend','painters','pos', [0 0 550 420]);
hold on
pressureProfile(tileA_all_0, stat, '*r');
pressureProfile(tileA_1_0, stat, '.k');
pressureProfile(tileA_2_0, stat, '.b');
pressureProfile(tileA_3_0, stat, '.g');
pressureProfile(tileA_4_0, stat, '.c');
pressureProfile(tileA_5_0, stat, '.m');
legend('all', '1', '2', '3', '4', '5');
pressureProfile(tileA_all_180, stat, '*r');
pressureProfile(tileA_1_180, stat, '.k');
pressureProfile(tileA_2_180, stat, '.b');
pressureProfile(tileA_3_180, stat, '.g');
pressureProfile(tileA_4_180, stat, '.c');
pressureProfile(tileA_5_180, stat, '.m');

return

% saveas(gcf, fullfile(path, sprintf('cp_%s_%s_00deg_profile_A.png', stat, s)))

% Tile B
figure('rend','painters','pos', [0 0 550 420]);
hold on
pressureProfile(tileB_1_0, stat, '*r');
pressureProfile(tileB_2_0, stat, 'ok');
legend('coarse', 'fine'); legend boxoff;
pressureProfile(tileB_1_180, stat, '*r');
pressureProfile(tileB_2_180, stat, 'ok');
% saveas(gcf, fullfile(path, sprintf('cp_%s_%s_00deg_profile_B.png', stat, s)))

return

% 20deg -------------------------------------------------------------------
figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_1_20, stat, ax);
hold on
pressureContour('no', tileB_1_20, stat, ax);
title('$$Cook - t = 18s, n = 16$$','interpreter','latex', 'fontsize', 22)
saveas(gcf, fullfile(path, sprintf('cppeak_%s_20deg_contour.png', 'time18')))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_2_20, stat, ax);
hold on
pressureContour('no', tileB_2_20, stat, ax);
title('$$Cook - t = 6s, n = 16$$','interpreter','latex', 'fontsize', 22)
saveas(gcf, fullfile(path, sprintf('cppeak_%s_20deg_contour.png', 'time6')))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_3_20, stat, ax);
hold on
pressureContour('no', tileB_3_20, stat, ax);
title('$$Cook - t = 2s, n = 16$$','interpreter','latex', 'fontsize', 22)
saveas(gcf, fullfile(path, sprintf('cppeak_%s_20deg_contour.png', 'time2')))

% Tile A
figure('rend','painters','pos', [0 0 550 420]);
hold on
pressureProfile(tileA_1_20, stat, '*r');
pressureProfile(tileA_2_20, stat, 'ok');
pressureProfile(tileA_3_20, stat, 'sb');
pressureProfile(tileA_0_20, stat, '--k');
legend('t = 18s, n = 16', 't = 6s, n = 16', 't = 2s, n = 16', 'safety factor'); legend boxoff;
saveas(gcf, fullfile(path, sprintf('cppeak_%s_20deg_profile_A.png', s)))

% Tile B
figure('rend','painters','pos', [0 0 550 420]);
hold on
pressureProfile(tileB_1_20, stat, '*r');
pressureProfile(tileB_2_20, stat, 'ok');
pressureProfile(tileB_3_20, stat, 'sb');
pressureProfile(tileB_0_20, stat, '--k');
legend('t = 18s, n = 16', 't = 6s, n = 16', 't = 2s, n = 16', 'safety factor'); legend boxoff;
saveas(gcf, fullfile(path, sprintf('cppeak_%s_20deg_profile_B.png', s)))


%% Pressure contours ------------------------------------------------------
% select plot
stat = 'peak'; % select variable to plot
ax = [-3 0];
path = '/home/giacomol/Desktop';

% 00deg -------------------------------------------------------------------
figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_poli_0, stat, ax);
hold on
pressureContour('no', tileB_poli_0, stat, ax);
pressureContour('no', tileA_poli_180, stat, ax);
pressureContour('no', tileB_poli_180, stat, ax);
title('$$PoliMi$$','interpreter','latex', 'fontsize', 22)
saveas(gcf, fullfile(path, 'cp_poli_peak_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_sens_0, stat, ax);
hold on
pressureContour('no', tileB_sens_0, stat, ax);
pressureContour('no', tileA_sens_180, stat, ax);
pressureContour('no', tileB_sens_180, stat, ax);
title('$$mean LES$$','interpreter','latex', 'fontsize', 22)
saveas(gcf, fullfile(path, 'cp_meanLES_peak_contour.png'))

return

ax = [0 1];
stat = 'spea';
figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_sens_0, stat, ax);
hold on
pressureContour('no', tileB_sens_0, stat, ax);
pressureContour('no', tileA_sens_180, stat, ax);
pressureContour('no', tileB_sens_180, stat, ax);
title('$$\Delta LES$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(path, 'cprms_LES_00deg_contour.png'))

return

% 20deg -------------------------------------------------------------------
figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_poli_20, stat, ax);
hold on
pressureContour('no', tileB_poli_20, stat, ax);
title('$$PoliMi$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(path, 'cprms_poli_20deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESk_20, stat, ax);
hold on
pressureContour('no', tileB_LESk_20, stat, ax);
title('$$LES$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(path, 'cprms_LES_20deg_contour.png'))

% 45deg -------------------------------------------------------------------
figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESk_45, stat, ax);
hold on
pressureContour('no', tileB_LESk_45, stat, ax);
title('$$LES$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(path, 'cprms_LES_45deg_contour.png'))

return

% RANS --------------------------------------------------------------------
figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', highRise_base_0, stat, ax);
title('$$C_{p} - central$$','interpreter','latex', 'fontsize', 22)
saveas(gcf, fullfile(path, 'cpmean_base_00deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', highRise_zhu_0, stat, ax);
title('$$C_{p}$$ - zhu','interpreter','latex', 'fontsize', 22)
saveas(gcf, fullfile(path, 'cpmean_zhu_00deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', highRise_upw_0, stat, ax);
title('$$C_{p} - upwind$$','interpreter','latex', 'fontsize', 22)
saveas(gcf, fullfile(path, 'cpmean_upwind_00deg_contour.png'))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', highRise_c2c_20, stat, ax);
title('$$C_{p}$$','interpreter','latex', 'fontsize', 22)

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', highRise_c3c_20, stat, ax);
title('$$C_{p}$$','interpreter','latex', 'fontsize', 22)

%% Pressure profiles ------------------------------------------------------
LES = [11,14,17];
stat = 'stan';

% find unique vertical coordinates on tile A
heights = unique(round(tileA_poli_0.coords(:,2),3));

% Tile A
figure('rend','painters','pos', [0 0 550 420]);
hold on
pressureProfile_all(highRise_base_0, 'A', stat, '-ok', 'k');
% pressureProfile_all(highRise_upw_0, 'A', stat, '-ob', 'b');
% pressureProfile_all(highRise_uref_0, 'A', stat, '-oc', 'c');
% pressureProfile_all(highRise_comb_20, 'A', stat, '.k', [0.5 0.5 0.5]);

% pressureProfile_all(highRise_LES_0, 'A', stat, '-sr', [0.5 0.5 0.5]); 

% pressureProfile_all(tileA_LESk{LES(1)}, 'A', stat, '-ob');
% pressureProfile_all(tileA_LESk{LES(2)}, 'A', stat, '-ok');
% pressureProfile_all(tileA_LESk{LES(3)}, 'A', stat, '-or');

% pressureProfile_all(tileA_sens, 'A', stat, '.k', [0.5 0.5 0.5]);

% legend({'$$z_0$$ - min','$$z_0$$ - mean','$$z_0$$ - max'}, 'interpreter', 'latex')
% legend({'$$k$$ - min','$$k$$ - mean','$$k$$ - max'}, 'interpreter', 'latex')
% legend({'$$T_u$$ - min','$$T_u$$ - mean','$$T_u$$ - max'}, 'interpreter', 'latex')

legend boxoff;

pressureProfile(tileA_poli_0, stat, 'or');
pressureProfile(tileA_poli_180, stat, 'or');
% ylim([-3 0])
% print(gcf, '-depsc','-painters', 'cp_epUQ_A.eps')
% epsclean('cp_epUQ_A.eps')
saveas(gcf, fullfile(path, sprintf('cp_A_peak_profile.png', stat)))

return

% Tile B
figure('rend','painters','pos', [0 0 550 420]);
hold on
% pressureProfile_all(highRise_base_0, 'B', stat, 'ok', 'k');
% pressureProfile_all(highRise_upw_0, 'B', stat, 'ob', 'b');
% pressureProfile_all(highRise_uref_0, 'B', stat, 'oc', 'c');
% pressureProfile_all(highRise_c3c_upw_20, 'B', stat, 'sc', 'c');
% pressureProfile_all(highRise_comb_20, 'B', stat, '.k', [0.5 0.5 0.5]);

pressureProfile(tileB_sens_0, stat, '.k', [0.5 0.5 0.5]); 
pressureProfile(tileB_sens_180, stat, '.k', [0.5 0.5 0.5]);

% pressureProfile_all(tileB_LESk_0{LES(1)}, 'B', stat, '-ob'); 
% pressureProfile_all(tileB_LESk_0{LES(2)}, 'B', stat, '-ok'); 
% pressureProfile_all(tileB_LESk_0{LES(3)}, 'B', stat, '-or');
pressureProfile(tileB_poli_0, stat, 'or');
pressureProfile(tileB_poli_180, stat, 'or');
% ylim([-3 0])

% legend('Smagorinsky', 'DynamicK', 'PoliMi'); legend boxoff
% print(gcf, '-depsc','-painters', 'cp_epUQ_B.eps')
% epsclean('cp_epUQ_B.eps')
saveas(gcf, fullfile(path, sprintf('cp_B_peak_profile.png', stat)))

return

%% accuracy
stat = 'peak';

[absA180, percA180] = accuracyCFD(tileA_sens_180, tileA_poli_180, stat);
[absA0, percA0] = accuracyCFD(tileA_sens_0, tileA_poli_0, stat);

[absB180, percB180] = accuracyCFD(tileB_sens_180, tileB_poli_180, stat);
[absB0, percB0] = accuracyCFD(tileB_sens_0, tileB_poli_0, stat);

err_abs = [absA0; absB0; absA180; absB180];
err_perc = [percA0; percB0; percA180; percB180];

Ntaps = 224*2 + 223*2;
sprintf('%6.4f%% of data captured', 100*(Ntaps - length(err_abs))/Ntaps)
sprintf('%6.4f%% mean percentage error', 100*mean(err_perc))
sprintf('%6.4f mean absolute error', mean(err_abs))

%% sensitivity profiles
figures = '/home/giacomol/Desktop/Research/slides/02-27-20/figures/'
params = {'z_0', 'tke', 'T_u'};
stat = 'desi'

for j = 1:3
    
    [tile_min, tile_mean, tile_max] = inflowAverageSensitivity(tileA_LES_180, stat, params{j});
    
    figure

    semilogx(tile_min.area, tile_min.design, '-ob')
    hold on
    semilogx(tile_mean.area, tile_mean.design, '-ok')
    semilogx(tile_max.area, tile_max.design, '-or')
    xlabel('$$A [m^2]$$', 'interpreter', 'latex')
    ylabel('$$C_{p,des}$$', 'interpreter', 'latex')
    set(gca,'fontsize',22)

    ylim([-1.5 -0.5])

    if params{j} == 'z_0'
        legend({'$$z_0$$ - min','$$z_0$$ - mean','$$z_0$$ - max'}, 'interpreter', 'latex', 'Location', 'NorthWest')
    elseif params{j}  == 'tke'
        legend({'$$k$$ - min','$$k$$ - mean','$$k$$ - max'}, 'interpreter', 'latex', 'Location', 'NorthWest')
    elseif params{j}  == 'T_u'
        legend({'$$T_u$$ - min','$$T_u$$ - mean','$$T_u$$ - max'}, 'interpreter', 'latex', 'Location', 'NorthWest')
    end
    legend boxoff;
    name = sprintf('cp_%s_%s_A.png', stat, params{j});
    ylim([-2.5 -0.5])
    saveas(gcf, fullfile(path, name))
end

return

%%
stats = {'peak'};
params =  {'z_0', 'tke', 'T_u'}
tile = 'B';

for i = 1:length(stats)
    for j = 1:length(params)
        if tile == 'A'
            [tile_min, tile_mean, tile_max] = inflowAverageSensitivity(tileA, stats{i}, params{j})

            figure('rend','painters','pos', [0 0 550 420]);
            hold on
            pressureProfile_all(tile_min,  tile, stats{i}, '-ob');
            pressureProfile_all(tile_mean, tile, stats{i}, '-ok');
            pressureProfile_all(tile_max,  tile, stats{i}, '-or');
            
        elseif tile == 'B'
            [tile_min_0,   tile_mean_0,   tile_max_0]   = inflowAverageSensitivity(tileB_LES_0, stats{i}, params{j})
            [tile_min_180, tile_mean_180, tile_max_180] = inflowAverageSensitivity(tileB_LES_180, stats{i}, params{j})
            
            figure('rend','painters','pos', [0 0 550 420]);
            hold on
            pressureProfile(tile_min_0,  stats{i}, '-ob');
            pressureProfile(tile_mean_0, stats{i}, '-ok');
            pressureProfile(tile_max_0,  stats{i}, '-or');     
            
            pressureProfile(tile_min_180,  stats{i}, '-ob');
            pressureProfile(tile_mean_180, stats{i}, '-ok');
            pressureProfile(tile_max_180,  stats{i}, '-or');             
        end

        if params{j} == 'z_0'
            legend({'$$z_0$$ - min','$$z_0$$ - mean','$$z_0$$ - max'}, 'interpreter', 'latex', 'Location', 'NorthWest')
        elseif params{j}  == 'tke'
            legend({'$$k$$ - min','$$k$$ - mean','$$k$$ - max'}, 'interpreter', 'latex', 'Location', 'NorthWest')
        elseif params{j}  == 'T_u'
            legend({'$$T_u$$ - min','$$T_u$$ - mean','$$T_u$$ - max'}, 'interpreter', 'latex', 'Location', 'NorthWest')
        end
        legend boxoff;
        name = sprintf('cp_%s_%s_%s.png', stats{i}, params{j}, tile);
        %ylim([0 0.3])
        saveas(gcf, fullfile(path, name))
    end
end

%% Sobol contours ---------------------------------------------------------
% LES sensitivity
[tileA_sens_0,   tileB_sens_0,   tileD_sens_0]   = pressureCombinedLES(tileA_0, tileB_0, tileD_LESk_0);
[tileA_sens_180, tileB_sens_180, tileD_sens_180] = pressureCombinedLES(tileA_180, tileB_180, tileD_LESk_180);

ax = [0 2];

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_sens_0, 'sob1', ax);
pressureContour('no', tileB_sens_0, 'sob1', ax);
pressureContour('no', tileA_sens_180, 'sob1', ax);
pressureContour('no', tileB_sens_180, 'sob1', ax);
title('$$S_{z_0}$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(figures, sprintf('s1_%s_meanUref.png', stat)))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_sens_0, 'sob2', ax);
pressureContour('no', tileB_sens_0, 'sob2', ax);
pressureContour('no', tileA_sens_180, 'sob2', ax);
pressureContour('no', tileB_sens_180, 'sob2', ax);
title('$$S_{k}$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(figures, sprintf('s2_%s_meanUref.png', stat)))

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_sens_0, 'sob3', ax);
pressureContour('no', tileB_sens_0, 'sob3', ax);
pressureContour('no', tileA_sens_180, 'sob3', ax);
pressureContour('no', tileB_sens_180, 'sob3', ax);
title('$$S_{T_u}$$','interpreter','latex', 'fontsize', 22)
% saveas(gcf, fullfile(figures, sprintf('s3_%s_meanUref.png', stat)))

return

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', highRise_base_0, 'sob1', [0 1]);
title('$$S_{z_0}$$','interpreter','latex', 'fontsize', 22)

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', highRise_base_0, 'sob2', [0 1]);
title('$$S_{U_{ref}}$$','interpreter','latex', 'fontsize', 22)

figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', highRise_base_0, 'sob3', [0 1]);
title('$$S_{\alpha}$$','interpreter','latex', 'fontsize', 22)

%% Sobol profiles ---------------------------------------------------------
% select plot
stat = 'mean'; % select variable to plot

% sobol
figure('rend','painters','pos', [0 0 550 420]);
pressureProfile(tileA_UQ_0, 'sob1', 'xk'); pressureProfile(tileA_UQ_180, 'sob1', 'xk');
hold on
pressureProfile(tileA_UQ_0, 'sob2', 'sr'); pressureProfile(tileA_UQ_180, 'sob2', 'sr');
pressureProfile(tileA_UQ_0, 'sob3', '.b'); pressureProfile(tileA_UQ_180, 'sob3', '.b');

% sobol
figure('rend','painters','pos', [0 0 550 420]);
pressureProfile(tileB_UQ_0, 'sob1', 'xk'); pressureProfile(tileB_UQ_180, 'sob1', 'xk');
hold on
pressureProfile(tileB_UQ_0, 'sob2', 'sr'); pressureProfile(tileB_UQ_180, 'sob2', 'sr');
pressureProfile(tileB_UQ_0, 'sob3', '.b'); pressureProfile(tileB_UQ_180, 'sob3', '.b');
        
%% Pressure time-series ---------------------------------------------------
list = 'A0301';

deltaT = 10;

% 00deg -------------------------------------------------------------------
figure('rend','painters','pos', [0 0 900 420]);

% PoliMi
tile = tileA_poli_180;
pressureTimeHistory('no', tile, 'r', -3, [1 tend], list);

% LES
tile = tileA_LESk_180{14};
tend = tile.time(end);
% [i,j] = find(tile.areaAverage == min(min(tile.areaAverage)));
% list = tile.taps(j,:)
% timeStep = tile.time(i);
pressureTimeHistory('no', tile, 'k', -3, [1 tend], list);
legend('PoliMi','LES')
legend boxoff

list = 'B0810';

figure('rend','painters','pos', [0 0 900 420]);

% PoliMi
tile = tileB_poli_180;
pressureTimeHistory('no', tile, 'r', -5, [1 tend], list);

% LES
tile = tileB_LESk_180{14};
pressureTimeHistory('no', tile, 'k', -5, [1 tend], list);

return

% 20deg -------------------------------------------------------------------
list = 'A0205';

% LES
tile = tileA_LESk_20;
tend = tile.time(end);
figure('rend','painters','pos', [0 0 900 420]);
pressureTimeHistory('no', tile, 'k', -5, [1 tend], list);

% PoliMi
tile = tileA_poli_20;
pressureTimeHistory('no', tile, 'r', -5, [1 tend], list);

list = 'B0810';

% LES
tile = tileB_LESk_20;
figure('rend','painters','pos', [0 0 900 420]);
pressureTimeHistory('no', tile, 'k', -5, [1 tend], list);

% PoliMi
tile = tileB_poli_20;
pressureTimeHistory('no', tile, 'r', -5, [1 tend], list);

%% Pressure spectra -------------------------------------------------------
list = 'A0205';

f = [0.001, 0.01, 0.1, 1, 100];

LES = [5,14,23];

figure
pressurePSD(tileA_poli_0, 500, 'adi', 'r', list);
hold on
pressurePSD(tileA_LES_0{14}, 1/0.0001, 'adi', 'k', list);
hold on
loglog(f, f.^(-2/3), '--k', 'linewidth', 1.5)
% legend('PoliMi','LES', 'Kolomogorov')
% legend boxoff

list = 'B0810';

figure
pressurePSD(tileB_poli_0, 500, 'adi', 'r', list);
hold on
pressurePSD(tileB_LES_0{14}, 1/0.0001, 'adi', 'k', list);
hold on
loglog(f, f.^(-2/3), '--k', 'linewidth', 1.5)

% legend({'$$PoliMi$$', '$$z_0$$ - min','$$z_0$$ - mean','$$z_0$$ - max'}, 'interpreter', 'latex')
% legend({'$$PoliMi$$', '$$k$$ - min','$$k$$ - mean','$$k$$ - max'}, 'interpreter', 'latex')
% legend({'$$PoliMi$$', '$$T_u$$ - min','$$T_u$$ - mean','$$T_u$$ - max'}, 'interpreter', 'latex')

%% pressure CDF -----------------------------------------------------------
list = 'B0808';

[~, ~, index] = intersect(list, highRise_base_20.taps(:,1:5), 'rows', 'stable');

% index = index + 224;

figure('rend','painters','pos', [0 0 550 420]);
hold on
plot(highRise_base_20.xcdf(:,index), highRise_base_20.fcdf(:,index), '--k', 'linewidth', 1.5)
plot(highRise_c1c_20.xcdf(:,index), highRise_c1c_20.fcdf(:,index))
plot(highRise_c2c_20.xcdf(:,index), highRise_c2c_20.fcdf(:,index))
plot(highRise_c3c_20.xcdf(:,index), highRise_c3c_20.fcdf(:,index))

[~, ~, index] = intersect(list, tileB_poli_20.taps(:,1:5), 'rows', 'stable');
plot(tileB_poli_20.mean(index)*ones(100,1), linspace(0,1), '--r')

ylabel('$$f$$','Interpreter','latex'); 
xlabel('$$C_P$$','Interpreter','latex');
set(gca,'fontsize', 18)
% axis([-0.5 0.0 0 1])

return

index = index + 224;

figure('rend','painters','pos', [0 0 550 420]);
hold on
plot(highRise_base_0.xcdf(:,index), highRise_base_0.fcdf(:,index), '--k', 'linewidth', 1.5)
plot(highRise_c1c_0.xcdf(:,index), highRise_c1c_0.fcdf(:,index))
plot(highRise_c2c_0.xcdf(:,index), highRise_c2c_0.fcdf(:,index))
plot(highRise_c3c_0.xcdf(:,index), highRise_c3c_0.fcdf(:,index))

[~, ~, index] = intersect(list, tileA_poli_180.taps(:,1:5), 'rows', 'stable');
plot(tileA_poli_180.mean(index)*ones(100,1), linspace(0,1), '--r')

ylabel('$$f$$','Interpreter','latex'); 
xlabel('$$C_P$$','Interpreter','latex');
set(gca,'fontsize', 18)
axis([-1 0.0 0 1])

figure
hold on
plotTAP3D(tileA_poli_0, '', 0)
% plotTAP3D(tileB_poli_0, '', 0)
%plotTAP3D(tileA_poli_180, list)
%plotTAP3D(tileB_poli_180, '')

%% Pressure PDF -----------------------------------------------------------
for i = 1:27

    figure(1)
    hold on
    pressurePDF(tileB_LESk_0{i}, [.5 .5 .5])

    figure(2)
    hold on
    pressurePDF(tileB_LESk_180{i}, [.5 .5 .5])   
end

figure(1)
pressurePDF(tileB_poli_0, [1 0 0])
% xlim([-2 1])

list = 'B0310';

figure(2)
pressurePDF(tileB_poli_180, [1 0 0])
% xlim([-2 1])

%% Extreme value analysis -------------------------------------------------
POE = 0.22;
window = 6;

for i = 1:27
    if i == 14
       col = [0 0 0]
    else
        col = [0.5 0.5 0.5]
    end
    figure(1); cp_peak_LES_A0(i) = gumbel(tileA_LES_0{i}, window, POE, col);
    figure(2); cp_peak_LES_B0(i) = gumbel(tileB_LES_0{i}, window, POE, col);
    figure(3); cp_peak_LES_A180(i) = gumbel(tileA_LES_180{i}, window, POE, col);
    figure(4); cp_peak_LES_B180(i) = gumbel(tileB_LES_180{i}, window, POE, col);
end
window = 6;
figure(1); cp_peak_poli_A0 = gumbel(tileA_poli_0, window, POE, [1 0 0]);
figure(2); cp_peak_poli_B0 = gumbel(tileB_poli_0, window, POE, [1 0 0]);
figure(3); cp_peak_poli_A180 = gumbel(tileA_poli_180, window, POE, [1 0 0]);
figure(4); cp_peak_poli_B180 = gumbel(tileB_poli_180, window, POE, [1 0 0]);


%% EVA vs area ------------------------------------------------------------
POE = 0.22;
ang = 180;

% PoliMi
[cp_peak_1, area] = designPressure('no', tileA_poli_180, ang, 1, POE, 'cook');
[cp_peak_2, area] = designPressure('no', tileA_poli_180, ang, 2, POE, 'cook');
[cp_peak_3, area] = designPressure('no', tileA_poli_180, ang, 3, POE, 'cook');
[cp_peak_4, area] = designPressure('no', tileA_poli_180, ang, 4, POE, 'cook');
[cp_peak_5, area] = designPressure('no', tileA_poli_180, ang, 5, POE, 'cook');

% WoW
% [cp_peak_wow, ~] = designPressure('no', tileA_wow_180, ang, 6, POE);

% LES
% for i = 1:27    
%     [cp_peak_LES(i,:), ~] = designPressure('no', tileA_LES_180{i}, ang, 6, POE, 'cook');
% end

%%
LES = [5,14,23];
figure
semilogx(area, cp_peak_5, '*r')
hold on
semilogx(area, cp_peak_1, '.')
semilogx(area, cp_peak_2, '.')
semilogx(area, cp_peak_3, '.')
semilogx(area, cp_peak_4, '.')
legend('all', '10 peaks', '20 peaks', '30 peaks', '40 peaks');
% semilogx(area, cp_peak_safe, 'ok')
% plot(area, cp_peak_wow, '.b', 'markersize', 15)
% errorbar(area, mean(cp_peak_LES), std(cp_peak_LES), 'ok')

% x = [area'; flipud(area')];
% y = [min(cp_peak_LES)'; flipud(max(cp_peak_LES)')];

% fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
 
% semilogx(area, cp_peak_LES(LES(1),:), '-ob')
% hold on
% semilogx(area, cp_peak_LES(LES(2),:), '-ok')
% semilogx(area, cp_peak_LES(LES(3),:), '-or')
set(gca,'fontsize',22)

% legend({'$$z_0$$ - min','$$z_0$$ - mean','$$z_0$$ - max'}, 'interpreter', 'latex')
% legend({'$$k$$ - min','$$k$$ - mean','$$k$$ - max'}, 'interpreter', 'latex')
% legend({'$$T_u$$ - min','$$T_u$$ - mean','$$T_u$$ - max'}, 'interpreter', 'latex')

legend boxoff;
set(gca, 'XScale', 'log');

xlabel('$$A [m^2]$$', 'interpreter','latex', 'fontsize', 22)
ylabel('$$C_{p,des}$$','interpreter','latex', 'fontsize', 22)
box off; %legend('PoliMi', 'WoW', 'LES'); legend boxoff
ylim([-2 -1])

return

S1 = sobol('z_0', cp_peak_LES');
S2 = sobol('TKE', cp_peak_LES');
S3 = sobol('T_u', cp_peak_LES');

figure
semilogx(area, S1, '-or')
hold on
semilogx(area, S2, '-sb')
semilogx(area, S3, '-xk')
box off; legend('S_{z_0}', 'S_{k}', 'S_{T_u}'); legend boxoff
xlabel('$$A [m^2]$$', 'interpreter','latex', 'fontsize', 22)
ylabel('$$S$$','interpreter','latex', 'fontsize', 22)
set(gca,'fontsize',18)
axis([5e-4 1e-1, 0 1.5])

%% Pressure instantaneous contours tileA ----------------------------------
% sharp peak --------------------------------------------------------------
T = 1000;
list = ['A0201';'A0502'];

% LES
tile  = tileA_LES_180{14};

[~, ~, index] = intersect(list, tile.taps(:,1:5), 'rows', 'stable');
tt_1 = find(tile.timeHistory(:,index(1)) <= -2.5);
t = tt_1(10);

figure
pressureInstantaneousContour(tile, tt_1(1), [-2.5, -1], list(1,:));
axis([0 0.1 1.9 2])
hold on


for i = 1:length(index)
    plot(tile.coords(index(i),1), tile.coords(index(i),2), 'ow', 'linewidth', 3)
end

figure
hold on
pressureTimeHistory('on', tile, '', -4, [tile.time(t-20*T), tile.time(t+20*T)], list(1,:));
pressureTimeHistory('on', tile, '', -5, [tile.time(t-20*T), tile.time(t+20*T)], list(2,:));
pressureTimeHistory('on', tile, '--k', -5, [tile.time(t-20*T), tile.time(t+20*T)]);
ylim([-3.5 0.5])

return

% PoliMi
tile = tileA_poli_180;

[~, ~, index] = intersect(list, tile.taps(:,1:5), 'rows', 'stable');
tt_1 = find(tile.timeHistory(:,index(1)) <= -2.75);
t = tt_1(1);

figure
pressureInstantaneousContour(tile, t, [-2.5, 0], list(1,:));
axis([0 0.1 1.9 2])
hold on

for i = 1:length(index)
    plot(tile.coords(index(i),1), tile.coords(index(i),2), 'ow', 'linewidth', 3)
end

figure
hold on
pressureTimeHistory('on', tile, '', -4, [tile.time(t-T), tile.time(t+T)], list(1,:));
pressureTimeHistory('on', tile, '', -5, [tile.time(t-T), tile.time(t+T)], list(2,:));
pressureTimeHistory('on', tile, '--k', -5, [tile.time(t-T), tile.time(t+T)]);
ylim([-3.5 0.5])

%% Pressure instantaneous contours tileB ----------------------------------
T = 1000;

% PoliMi
list = ['B1203';'B0310'];
tile = tileB_poli_180;

[~, ~, index] = intersect(list, tile.taps(:,1:5), 'rows', 'stable');
tt_1 = find(tile.timeHistory(:,index(1)) <= -2);
t = tt_1(1);

figure
pressureInstantaneousContour(tile, t, [-2.5, -1], list(1,:));
% axis([0 0.1 1.9 2])
hold on

for i = 1:length(index)
    plot(tile.coords(index(i),1), tile.coords(index(i),2), 'ow', 'linewidth', 3)
end
% 
% figure
% hold on
% pressureTimeHistory('on', tile, '', -4, [tile.time(t-T), tile.time(t+T)], list(1,:));
% pressureTimeHistory('on', tile, '', -5, [tile.time(t-T), tile.time(t+T)], list(2,:));
% pressureTimeHistory('on', tile, '--k', -5, [tile.time(t-T), tile.time(t+T)]);
% ylim([-3 0])

% LES
tile  = tileB_LESk_180{14};

[~, ~, index] = intersect(list, tile.taps(:,1:5), 'rows', 'stable');
tt_1 = find(tile.timeHistory(:,index(1)) <= -1.5);
t = tt_1(end);

figure
pressureInstantaneousContour(tile, tt_1(1), [-2.5, -1], list(1,:));
% axis([0 0.1 1.9 2])
hold on


for i = 1:length(index)
    plot(tile.coords(index(i),1), tile.coords(index(i),2), 'ow', 'linewidth', 3)
end
% 
% figure
% hold on
% pressureTimeHistory('on', tile, '', -4, [tile.time(t-20*T), tile.time(t+20*T)], list(1,:));
% pressureTimeHistory('on', tile, '', -5, [tile.time(t-20*T), tile.time(t+20*T)], list(2,:));
% pressureTimeHistory('on', tile, '--k', -5, [tile.time(t-20*T), tile.time(t+20*T)]);
% ylim([-3 0])

%% movie ------------------------------------------------------------------
% PoliMi
writerObj = VideoWriter('cp_contour+timeHistory_poli_A180.avi');
writerObj.FrameRate = 4;
tile  = tileA_poli_180;
[i,j] = find(tile.timeHistory == min(min(tile.timeHistory)));
list  = [tile.taps(j,:); 'A0808'];
pressureMovie(tile, i, 0.01, writerObj, [-2 0], list)

% LES
writerObj = VideoWriter('cp_contour+timeHistory_smagorinsky_A180.avi');
writerObj.FrameRate = 4;
tile  = tileA_LES_180{14};
[i,j] = find(tile.timeHistory == min(min(tile.timeHistory)));
list  = [tile.taps(j,:); 'A0808'];
pressureMovie(tile, i, 0.01, writerObj, [-2 0], list)

% LES
writerObj = VideoWriter('cp_contour+timeHistory_dynamicK_A180.avi');
writerObj.FrameRate = 4;
tile  = tileA_LESk_180{14};
[i,j] = find(tile.timeHistory == min(min(tile.timeHistory)));
list  = [tile.taps(j,:); 'A0808'];
pressureMovie(tile, i, 0.01, writerObj, [-2 0], list)

%% SPOD -------------------------------------------------------------------
[phi, lambda, Pxx] = computeSPODmodes(tileA_LESk_0{14}, 0.0001, 10, 100);
tileA_LESk_0{14}.spod = phi;

[phi, lambda, Pxx] = computeSPODmodes(tileB_LESk_0{14}, 0.0001, 10, 100);
tileB_LESk_0{14}.spod = phi;

[phi, lambda, Pxx] = computeSPODmodes(tileA_LESk_180{14}, 0.0001, 10, 100);
tileA_LESk_180{14}.spod = phi;

[phi, lambda, Pxx] = computeSPODmodes(tileB_LESk_180{14}, 0.0001, 10, 100);
tileB_LESk_180{14}.spod = phi;

%%
figure('rend','painters','pos', [0 0 400 400]);
pressureContour('no', tileA_LESk_0{14}, 'spod', ax);
pressureContour('no', tileB_LESk_0{14}, 'spod', ax);
pressureContour('no', tileA_LESk_180{14}, 'spod', ax);
pressureContour('no', tileB_LESk_180{14}, 'spod', ax);
title('$$S_{z_0}$$','interpreter','latex', 'fontsize', 22)
