%% Compare CharLES and Wind tunnel design pressure

clear; close all; clc

angle = 10;
Tcharles = 98;
Texp = 300;

%% Load charLES data
tileA_charles = load_charLES(angle, 'A', sprintf('/home/giacomol/Desktop/Research/windLoading/charLES_fine/%d/pcprobespA0/', angle));
tileB_charles = load_charLES(angle, 'B', sprintf('/home/giacomol/Desktop/Research/windLoading/charLES_fine/%d/pcprobespB0/', angle));

%% Load wind tunnel data
% PoliMi
[tileA_poli, tileB_poli, ~] = pressurePoliMi(angle);

% WoW
[tileA_wow, tileB_wow, ~] = pressureWoW(angle);

%% Design pressure
dx = 0.06;
dy = 0.04; 

window = 6;
start = 0;

path = '/home/giacomol/Desktop/Research/slides/10-29-20/figures/';

% 20deg
figure('Position', [10 10 900 600])
tileA_poli = designPressurePanels(tileA_poli, dx, dy, window, floor(tileA_poli.time(end)/window), start);
tileB_poli = designPressurePanels(tileB_poli, dx, dy, window, floor(tileB_poli.time(end)/window), start);
title('PoliMi')
saveas(gcf, fullfile(path, sprintf('cp_design_poli_%ddeg_window6s.png', angle)))

figure('Position', [10 10 900 600])
tileA_wow = designPressurePanels(tileA_wow, dx, dy, window, floor(tileA_wow.time(end)/window), start); 
tileB_wow = designPressurePanels(tileB_wow, dx, dy, window, floor(tileB_wow.time(end)/window), start); 
title('WoW')
saveas(gcf, fullfile(path, sprintf('cp_design_wow_%ddeg_window6s.png', angle)))

figure('Position', [10 10 900 600])
tileA_charles = designPressurePanels(tileA_charles, dx, dy, window, floor(tileA_charles.time(end)/window), start); 
tileB_charles = designPressurePanels(tileB_charles, dx, dy, window, floor(tileA_charles.time(end)/window), start); 
title('CharLES')
saveas(gcf, fullfile(path, sprintf('cp_design_charles_%ddeg_window6s.png', angle)))

% 20deg
figure
plot(tileA_charles.design, tileA_poli.design, '.r', 'markersize', 10)
hold on
plot(tileA_charles.design, tileA_wow.design, '.b', 'markersize', 10)
plot(linspace(-3,0), linspace(-3,0), '--k')
xlabel('$$C_{p,des,LES}$$','interpreter','latex')
ylabel('$$C_{p,des,WT}$$','interpreter','latex')
legend({'PoliMi', 'WoW'}, 'location', 'southeast')
set(gca, 'fontsize', 18)
axis([-2.5, -0.5, -2.5, -0.5])
pbaspect([1 1 1])
saveas(gcf, fullfile(path, sprintf('fit_A%ddeg_window6s.png', angle)))

figure
plot(tileB_charles.design, tileB_poli.design, '.r', 'markersize', 10)
hold on
plot(tileB_charles.design, tileB_wow.design, '.b', 'markersize', 10)
plot(linspace(-3,0), linspace(-3,0), '--k')
xlabel('$$C_{p,des,LES}$$','interpreter','latex')
ylabel('$$C_{p,des,WT}$$','interpreter','latex')
legend({'PoliMi', 'WoW'}, 'location', 'southeast')
set(gca, 'fontsize', 18)
axis([-2.5, -0.5, -2.5, -0.5])
pbaspect([1 1 1])
saveas(gcf, fullfile(path, sprintf('fit_B%ddeg_window6s.png', angle)))

%% Fix number of peaks
npeaks = 16;

k = 0;
for start = 0 : floor(tileA_poli.time(end)/window) - npeaks
    
    k = k + 1;    
    tile = designPressurePanels(tileA_poli, dx, dy, window, npeaks, start);
    cp_design_poli_A(:,k) = tile.design;
    
    tile = designPressurePanels(tileB_poli, dx, dy, window, npeaks, start);
    cp_design_poli_B(:,k) = tile.design;
end

k = 0;
for start = 0 : floor(tileA_wow.time(end)/window) - npeaks
    
    k = k + 1;    
    tile = designPressurePanels(tileA_wow, dx, dy, window, npeaks, start); 
    cp_design_wow_A(:,k) = tile.design;
    
    tile = designPressurePanels(tileB_wow, dx, dy, window, npeaks, start);
    cp_design_wow_B(:,k) = tile.design;
end

%%
% 20deg
figure
errorbar(tileA_charles.design, mean(cp_design_poli_A'), 1.96 * std(cp_design_poli_A'), 'or')
hold on
errorbar(tileA_charles.design, mean(cp_design_wow_A'), 1.96 * std(cp_design_wow_A'), 'ob')
plot(linspace(-3,0), linspace(-3,0), '--k')
xlabel('$$C_{p,des,LES}$$','interpreter','latex')
ylabel('$$C_{p,des,WT}$$','interpreter','latex')
legend({'PoliMi', 'WoW'}, 'location', 'southeast')
set(gca, 'fontsize', 18)
axis([-2.5, -0.5, -2.5, -0.5])
pbaspect([1 1 1])
saveas(gcf, fullfile(path, sprintf('fit_A%ddeg_npeaks16_vs_windows.png', angle)))

figure
errorbar(tileB_charles.design, mean(cp_design_poli_B'), 1.96 * std(cp_design_poli_B'), 'or')
hold on
errorbar(tileB_charles.design, mean(cp_design_wow_B'), 1.96 * std(cp_design_wow_B'), 'ob')
plot(linspace(-3,0), linspace(-3,0), '--k')
xlabel('$$C_{p,des,LES}$$','interpreter','latex')
ylabel('$$C_{p,des,WT}$$','interpreter','latex')
legend({'PoliMi', 'WoW'}, 'location', 'southeast')
set(gca, 'fontsize', 18)
axis([-2.5, -0.5, -2.5, -0.5])
pbaspect([1 1 1])
saveas(gcf, fullfile(path, sprintf('fit_B%ddeg_npeaks16_vs_windows.png', angle)))

clear cp_design_poli_A cp_design_poli_B
clear cp_design_wow_A cp_design_wow_B
clear cp_design_charles_A cp_design_charles_B

%%
npeaks = 16;
start = 0;
k = 0;

for window = 1 : floor(Texp/npeaks)
    
    k = k + 1;
    tile = designPressurePanels(tileA_poli, dx, dy, window, npeaks, start);
    cp_design_poli_A(:,k) = tile.design;
    
    tile = designPressurePanels(tileB_poli, dx, dy, window, npeaks, start);
    cp_design_poli_B(:,k) = tile.design;
    
    tile = designPressurePanels(tileA_wow, dx, dy, window, npeaks, start); 
    cp_design_wow_A(:,k) = tile.design;
    
    tile = designPressurePanels(tileB_wow, dx, dy, window, npeaks, start);
    cp_design_wow_B(:,k) = tile.design;
    
    if window <= floor(Tcharles/npeaks)
        tile = designPressurePanels(tileA_charles, dx, dy, window, npeaks, start);
        cp_design_charles_A(:,k) = tile.design;

        tile = designPressurePanels(tileB_charles, dx, dy, window, npeaks, start);
        cp_design_charles_B(:,k) = tile.design;     
    end
end

figure
hold on
plot(1 : floor(Texp/npeaks), cp_design_poli_A(150,:), '-ob')
plot(1 : floor(Texp/npeaks), cp_design_wow_A(150,:), '-or')
plot(1 : floor(Tcharles/npeaks), cp_design_charles_A(150,:), '-ok')
legend({'PoliMi', 'WoW', 'CharLES'}, 'location', 'northeast')
xlabel('$$Window$$','interpreter','latex')
ylabel('$$C_{p,des}$$','interpreter','latex')
set(gca, 'fontsize', 18)
ylim([-4 -0.5])
saveas(gcf, fullfile(path, sprintf('cp_design_A%ddeg_panel150_npeaks16_vs_windows.png', angle)))

figure
hold on
plot(1 : floor(Texp/npeaks), cp_design_poli_A(65,:), '-ob')
plot(1 : floor(Texp/npeaks), cp_design_wow_A(65,:), '-or')
plot(1 : floor(Tcharles/npeaks), cp_design_charles_A(65,:), '-ok')
legend({'PoliMi', 'WoW', 'CharLES'}, 'location', 'southeast')
xlabel('$$Window$$','interpreter','latex')
ylabel('$$C_{p,des}$$','interpreter','latex')
set(gca, 'fontsize', 18)
ylim([-3 -0.5])
saveas(gcf, fullfile(path, sprintf('cp_design_A%ddeg_panel65_npeaks16_vs_windows.png', angle)))

figure
hold on
plot(1 : floor(Texp/npeaks), cp_design_poli_B(40,:), '-ob')
plot(1 : floor(Texp/npeaks), cp_design_wow_B(40,:), '-or')
plot(1 : floor(Tcharles/npeaks), cp_design_charles_B(40,:), '-ok')
legend({'PoliMi', 'WoW', 'CharLES'}, 'location', 'northeast')
xlabel('$$Window$$','interpreter','latex')
ylabel('$$C_{p,des}$$','interpreter','latex')
set(gca, 'fontsize', 18)
ylim([-3 -0.5])
saveas(gcf, fullfile(path, sprintf('cp_design_B%ddeg_panel40_npeaks16_vs_windows.png', angle)))

figure
hold on
plot(1 : floor(Texp/npeaks), cp_design_poli_B(17,:), '-ob')
plot(1 : floor(Texp/npeaks), cp_design_wow_B(17,:), '-or')
plot(1 : floor(Tcharles/npeaks), cp_design_charles_B(17,:), '-ok')
legend({'PoliMi', 'WoW', 'CharLES'}, 'location', 'southeast')
xlabel('$$Window$$','interpreter','latex')
ylabel('$$C_{p,des}$$','interpreter','latex')
set(gca, 'fontsize', 18)
ylim([-3 -0.5])
saveas(gcf, fullfile(path, sprintf('cp_design_B%ddeg_panel17_npeaks16_vs_windows.png', angle)))

clear cp_design_poli_A cp_design_poli_B
clear cp_design_wow_A cp_design_wow_B
clear cp_design_charles_A cp_design_charles_B

%%
Texp = 300;
window = 6;
start = 0;
k = 0;

for npeaks = 10 : floor(Texp/window) - 10
    
    k = k + 1;
    tile = designPressurePanels(tileA_poli, dx, dy, window, npeaks, start);
    cp_design_poli_A(:,k) = tile.design;
    
    tile = designPressurePanels(tileB_poli, dx, dy, window, npeaks, start);
    cp_design_poli_B(:,k) = tile.design;
    
    tile = designPressurePanels(tileA_wow, dx, dy, window, npeaks, start); 
    cp_design_wow_A(:,k) = tile.design;
    
    tile = designPressurePanels(tileB_wow, dx, dy, window, npeaks, start);
    cp_design_wow_B(:,k) = tile.design;
    
    if npeaks <= floor(Tcharles/window)
        tile = designPressurePanels(tileA_charles, dx, dy, window, npeaks, start);
        cp_design_charles_A(:,k) = tile.design;

        tile = designPressurePanels(tileB_charles, dx, dy, window, npeaks, start);
        cp_design_charles_B(:,k) = tile.design;     
    end    
end

figure
hold on
plot(10 : floor(Texp/window) - 10, cp_design_poli_A(150,:), '-ob')
plot(10 : floor(Texp/window) - 10, cp_design_wow_A(150,:), '-or')
plot(10 : floor(Tcharles/window), cp_design_charles_A(150,:), '-ok')
legend({'PoliMi', 'WoW', 'CharLES'}, 'location', 'northeast')
xlabel('$$N peaks$$','interpreter','latex')
ylabel('$$C_{p,des}$$','interpreter','latex')
set(gca, 'fontsize', 18)
ylim([-4 -0.5])
saveas(gcf, fullfile(path, sprintf('cp_design_A%ddeg_panel150_npeaks_vs_window6.png', angle)))

figure
hold on
plot(10 : floor(Texp/window) - 10, cp_design_poli_A(65,:), '-ob')
plot(10 : floor(Texp/window) - 10, cp_design_wow_A(65,:), '-or')
plot(10 : floor(Tcharles/window), cp_design_charles_A(65,:), '-ok')
legend({'PoliMi', 'WoW', 'CharLES'}, 'location', 'southeast')
xlabel('$$N peaks$$','interpreter','latex')
ylabel('$$C_{p,des}$$','interpreter','latex')
set(gca, 'fontsize', 18)
ylim([-3 -0.5])
saveas(gcf, fullfile(path, sprintf('cp_design_A%ddeg_panel65_npeaks_vs_window6.png', angle)))

figure
hold on
plot(10 : floor(Texp/window) - 10, cp_design_poli_B(40,:), '-ob')
plot(10 : floor(Texp/window) - 10, cp_design_wow_B(40,:), '-or')
plot(10 : floor(Tcharles/window), cp_design_charles_B(40,:), '-ok')
legend({'PoliMi', 'WoW', 'CharLES'}, 'location', 'northeast')
xlabel('$$N peaks$$','interpreter','latex')
ylabel('$$C_{p,des}$$','interpreter','latex')
set(gca, 'fontsize', 18)
ylim([-3 -0.5])
saveas(gcf, fullfile(path, sprintf('cp_design_B%ddeg_panel40_npeaks_vs_window6.png', angle)))

figure
hold on
plot(10 : floor(Texp/window) - 10, cp_design_poli_B(17,:), '-ob')
plot(10 : floor(Texp/window) - 10, cp_design_wow_B(17,:), '-or')
plot(10 : floor(Tcharles/window), cp_design_charles_B(17,:), '-ok')
legend({'PoliMi', 'WoW', 'CharLES'}, 'location', 'southeast')
xlabel('$$N peaks$$','interpreter','latex')
ylabel('$$C_{p,des}$$','interpreter','latex')
set(gca, 'fontsize', 18)
ylim([-3 -0.5])
saveas(gcf, fullfile(path, sprintf('cp_design_B%ddeg_panel17_npeaks_vs_window6.png', angle)))
