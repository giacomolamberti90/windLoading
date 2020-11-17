%% POST PROCESSING LES

clear; close all; clc

% case directory and file
path_EXP = '/home/giacomol/Desktop/Research/windLoading/windTunnel/pythonCode/exp';

% velocity matrices
z_exp = dlmread(fullfile(path_EXP, 'UInlet'), '', [1  0  49  0]);

U_exp_mean = dlmread(fullfile(path_EXP, 'UInlet'), '', [1  1  49  1]);
U_exp_min  = dlmread(fullfile(path_EXP, 'UInlet_min'), '', [1  1  49  1]);
U_exp_max  = dlmread(fullfile(path_EXP, 'UInlet_max'), '', [1  1  49  1]);

uu_exp_mean = dlmread(fullfile(path_EXP, 'uuBarInlet'), '', [1  1  49  1]);
uu_exp_min  = dlmread(fullfile(path_EXP, 'uuBarInlet_min'), '', [1  1  49  1]);
uu_exp_max  = dlmread(fullfile(path_EXP, 'uuBarInlet_max'), '', [1  1  49  1]);

vv_exp_mean = dlmread(fullfile(path_EXP, 'vvBarInlet'), '', [1  1  49  1]);
vv_exp_min  = dlmread(fullfile(path_EXP, 'vvBarInlet_min'), '', [1  1  49  1]);
vv_exp_max  = dlmread(fullfile(path_EXP, 'vvBarInlet_max'), '', [1  1  49  1]);

ww_exp_mean = dlmread(fullfile(path_EXP, 'wwBarInlet'), '', [1  1  49  1]);
ww_exp_min  = dlmread(fullfile(path_EXP, 'wwBarInlet_min'), '', [1  1  49  1]);
ww_exp_max  = dlmread(fullfile(path_EXP, 'wwBarInlet_max'), '', [1  1  49  1]);

xLu_exp_mean = dlmread(fullfile(path_EXP, 'xLuInlet'), '', [1  1  49  1]);
xLu_exp_min  = dlmread(fullfile(path_EXP, 'xLuInlet_min'), '', [1  1  49  1]);
xLu_exp_max  = dlmread(fullfile(path_EXP, 'xLuInlet_max'), '', [1  1  49  1]);

Tu_exp_mean = xLu_exp_mean./U_exp_mean;
Tu_exp_min  = xLu_exp_min./U_exp_max;
Tu_exp_max  = xLu_exp_max./U_exp_min;

yLu_exp_mean = dlmread(fullfile(path_EXP, 'yLuInlet'), '', [1  1  49  1]);
yLu_exp_min  = dlmread(fullfile(path_EXP, 'yLuInlet_min'), '', [1  1  49  1]);
yLu_exp_max  = dlmread(fullfile(path_EXP, 'yLuInlet_max'), '', [1  1  49  1]);

zLu_exp_mean = dlmread(fullfile(path_EXP, 'zLuInlet'), '', [1  1  49  1]);
zLu_exp_min  = dlmread(fullfile(path_EXP, 'zLuInlet_min'), '', [1  1  49  1]);
zLu_exp_max  = dlmread(fullfile(path_EXP, 'zLuInlet_max'), '', [1  1  49  1]);

% TKE
k_exp_mean = 0.5 * (uu_exp_mean + vv_exp_mean + ww_exp_mean);
k_exp_min  = 0.5 * (uu_exp_min + vv_exp_min + ww_exp_min);
k_exp_max  = 0.5 * (uu_exp_max + vv_exp_max + ww_exp_max);

% turbulence intensities
Iu_exp_mean = sqrt(uu_exp_mean)./U_exp_mean;
Iu_exp_min  = sqrt(uu_exp_min)./U_exp_max;
Iu_exp_max  = sqrt(uu_exp_max)./U_exp_min;

Iv_exp_mean = sqrt(vv_exp_mean)./U_exp_mean;
Iv_exp_min  = sqrt(vv_exp_min)./U_exp_max;
Iv_exp_max  = sqrt(vv_exp_max)./U_exp_min;

Iw_exp_mean = sqrt(ww_exp_mean)./U_exp_mean;
Iw_exp_min  = sqrt(ww_exp_min)./U_exp_max;
Iw_exp_max  = sqrt(ww_exp_max)./U_exp_min;

%% LES --------------------------------------------------------------------
path_SENS = '/home/giacomol/Desktop/Research/windLoading/LES/highRise/00deg/coarse';
Nprobes = 21;
NLES = 3;
Nback = 5;
Nwindows = 1;

z = linspace(0, 4, Nprobes);

for num = 1%:NLES
    
    if num == 1
        path_LES = fullfile(path_SENS, 'dynamicK/workdir.14');
    elseif num == 2
        path_LES = fullfile(path_SENS, 'dynamicLagrangian/workdir.14');
    elseif num == 3
        path_LES = fullfile(path_SENS, 'smagorinsky/workdir.14');
    end
    
    U{num}  = dlmread(fullfile(path_LES, '/U_x_2m.out'));    
    
    uu{num} = dlmread(fullfile(path_LES, '/uu_z1p75m.out'));
    vv{num} = dlmread(fullfile(path_LES, '/vv_z1p75m.out'));
    ww{num} = dlmread(fullfile(path_LES, '/ww_z1p75m.out'));
    
    Iu{num} = sqrt(uu{num})./U{num};
    Iv{num} = sqrt(vv{num})./U{num};
    Iw{num} = sqrt(ww{num})./U{num};
    
    TKE{num} = 0.5 * (uu{num} + vv{num} + ww{num});
    
    % time-scales
    % time in probes
    dire = dir(fullfile(path_LES, sprintf('/postProcessing/probes_z1p75/')));
    for i = 1:length(dire)-2
        tprobes(i) = str2num(dire(2+i).name);
    end
    tprobes = sort(tprobes);
    
    data_L_1p75 = [];
    Nback = length(tprobes)-1;
    for i = length(tprobes)-Nback:length(tprobes)
        % probes -> length-scales
        if num == 14
            if floor(tprobes(i)) == tprobes(i)
                fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2d/U', tprobes(i))));
            else
                fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2.2f/U', tprobes(i))));
            end            
        else
            if floor(tprobes(i)) == tprobes(i)
                fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2d/U', tprobes(i))));
            else
                fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2.1f/U', tprobes(i))));
            end
        end
        format  = repmat('%f', 1, (Nprobes+1)*3);
        data_Li_1p75 = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes+2);
        fclose(fid);
        
        data_Li_1p75 = cell2mat(data_Li_1p75);
        
        if i == length(tprobes)-Nback
            tin = 1;
        else
            tin = find(data_Li_1p75(:,1) == tend) + 1;
        end
        data_L_1p75 = [data_L_1p75; data_Li_1p75(tin:end,:)];
        tend = data_Li_1p75(end,1);
    end
    
    % time
    time = data_L_1p75(:,1);
    N    = length(time);
        
    for i = 1:Nprobes
        u_1p75(:,i) = data_L_1p75(:,2+3*(i-1));
        Tu_1p75{num}(i) = computeTimeScale(time-time(1), data_L_1p75(:,2+3*(i-1)), Nwindows);
        %Lu_1p75{num}(i) = Tu_1p75{num}(i);
    end
    clear data_Li_1p75 data_L_1p75 tprobes
    
    % time-scales
    dire = dir(fullfile(path_LES, sprintf('/postProcessing/probes_z_1p45/')));
    for i = 1:length(dire)-2
        tprobes(i) = str2num(dire(2+i).name);
    end
    tprobes = sort(tprobes);
    
    data_L_1p45 = [];
    Nback = length(tprobes)-1;
    for i = length(tprobes)-Nback:length(tprobes)
        % probes -> length-scales
        if num == 14
            if floor(tprobes(i)) == tprobes(i)
                fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2d/U', tprobes(i))));
            else
                fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2.2f/U', tprobes(i))));
            end            
        else
            if floor(tprobes(i)) == tprobes(i)
                fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2d/U', tprobes(i))));
            else
                fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2.1f/U', tprobes(i))));
            end
        end
        format  = repmat('%f', 1, (Nprobes+1)*3);
        data_Li_1p45 = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes+2);
        fclose(fid);
        
        data_Li_1p45 = cell2mat(data_Li_1p45);
        
        if i == length(tprobes)-Nback
            tin = 1;
        else
            tin = find(data_Li_1p45(:,1) == tend) + 1;
        end
        data_L_1p45 = [data_L_1p45; data_Li_1p45(tin:end,:)];
        tend = data_Li_1p45(end,1);
    end
    
    % time
    time = data_L_1p45(:,1);
    N    = length(time);
        
    for i = 1:Nprobes
        u_1p45(:,i) = data_L_1p45(:,2+3*(i-1));
        Tu_1p45{num}(i) = computeTimeScale(time-time(1), data_L_1p45(:,2+3*(i-1)), Nwindows);
        %Lu_1p45{num}(i) = Tu_1p45{num}(i);
    end
    clear data_Li_1p45 data_L_1p45 tprobes
    
    Tu{num} = (Tu_1p45{num} + Tu_1p75{num})/2;
    %Lu{num} = (Lu_1p45{num} + Lu_1p75{num})/2;
end

%% figures
col = {'k','c','b'};
path = '/home/giacomol/Desktop/';

path_RANS = '/home/giacomol/Desktop/Research/windLoading/RANS/highRise/ML/00deg';
U_RANS = dlmread(fullfile(path_RANS, '/U_x_4m.out'));
k_RANS = dlmread(fullfile(path_RANS, '/k_x_4m.out'));

% U
figure
hold on
set(gca, 'fontsize', 18)
for num = 1%:NLES
    plot(U{num}, z, '-ok', 'linewidth', 2)
end
plot(U_RANS, z, '--b', 'linewidth', 2)
% plot(U_exp_mean, z_exp, '--k')
% % legend('DynamicK', 'Lagrangian', 'Smagorinsky', 'PoliMi'); legend boxoff
% y = [z_exp; flipud(z_exp)];
% x = [U_exp_min; flipud(U_exp_max)];
% fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$U [m/s]$$', 'interpreter', 'latex');
axis([0 9 0 4])
% saveas(gcf, fullfile(path, 'U_LES+RANS.png'))

% R11
figure
hold on
set(gca, 'fontsize', 18)
for num = 1%:NLES
    plot(Iu{num}, z, '-ok', 'linewidth', 2)
end
% plot(Iu_exp_mean, z_exp, '--k')
% % legend('DynamicK', 'Lagrangian', 'Smagorinsky', 'PoliMi'); legend boxoff
% y = [z_exp; flipud(z_exp)];
% x = [Iu_exp_min; flipud(Iu_exp_max)];
% fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$I_u$$', 'interpreter', 'latex');
axis([0 0.4 0 4])
% saveas(gcf, fullfile(path, 'Iu_LES+RANS.png'))

% R22
figure
hold on
set(gca, 'fontsize', 18)
for num = 1%:NLES
    plot(Iv{num}, z, '-ok', 'linewidth', 2)
end
% plot(Iv_exp_mean, z_exp, '--k')
% % legend('DynamicK', 'Lagrangian', 'Smagorinsky', 'PoliMi'); legend boxoff
% y = [z_exp; flipud(z_exp)];
% x = [Iv_exp_min; flipud(Iv_exp_max)];
% fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$I_v$$', 'interpreter', 'latex');
axis([0 0.3 0 4])
% saveas(gcf, fullfile(path, 'Iv_LES+RANS.png'))

% R33
figure
hold on
set(gca, 'fontsize', 18)
for num = 1%:NLES
    plot(Iw{num}, z, '-ok', 'linewidth', 2)
end
% plot(Iw_exp_mean, z_exp, '--k')
% % legend('DynamicK', 'Lagrangian', 'Smagorinsky', 'PoliMi'); legend boxoff
% y = [z_exp; flipud(z_exp)];
% x = [Iw_exp_min; flipud(Iw_exp_max)];
% fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$I_w$$', 'interpreter', 'latex');
axis([0 0.3 0 4])
% saveas(gcf, fullfile(path, 'Iw_LES+RANS.png'))

% TKE
figure
hold on
set(gca, 'fontsize', 18)
for num = 1%:NLES
    plot(TKE{num}, z, '-ok', 'linewidth', 2)
end
plot(k_RANS, z, '--b', 'linewidth', 2)
% plot(k_exp_mean, z_exp, '--k')
% % legend('DynamicK', 'Lagrangian', 'Smagorinsky', 'PoliMi'); legend boxoff
% y = [z_exp; flipud(z_exp)];
% x = [k_exp_min; flipud(k_exp_max)];
% fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$k [m^2/s^2]$$', 'interpreter', 'latex');
axis([0 2 0 4])
% saveas(gcf, fullfile(path, 'k_LES+RANS.png'))

% Tu
figure
hold on
set(gca, 'fontsize', 18)
for num = 1%:NLES
    plot(Tu{num}, z, '-ok', 'linewidth', 2)
end
% plot(Tu_exp_mean, z_exp, '--k')
% % legend('DynamicK', 'Lagrangian', 'Smagorinsky', 'PoliMi'); legend boxoff
% y = [z_exp; flipud(z_exp)];
% x = [Tu_exp_min; flipud(Tu_exp_max)];
% fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$T_u [s]$$', 'interpreter', 'latex');
axis([0 0.1 0 4])
% saveas(gcf, fullfile(path, 'Tu_LES+RANS.png'))

%%
idx = 11;

fsamp = 1/(time(2) - time(1));
[S, f] = pwelch(u_1p75(:,idx) - mean(u_1p75(:,idx)), [], [], [], fsamp);

K53 = f(2:end).^(-2/3);

figure
loglog(f(2:end), S(2:end), 'linewidth', 1.5, 'color', 'k')
hold on
loglog(f(2:end), 1.5 * K53, 'r', 'linewidth', 1.5)

Lu = mean(Tu{1}(idx)) * U{1}(idx);
loglog(f(2:end)*Lu/U{1}(idx), 4*(f(2:end)*Lu/U{1}(idx))./(1+70.8*(f(2:end)*Lu/U{1}(idx)).^2).^(5/6), 'b')

xlabel('$$f [Hz]$$','Interpreter','latex'); 
ylabel('$$S [m^2/s^2]$$','Interpreter','latex');

axis([0.01, 1000, 1e-6, 10]);

set(gca, 'fontsize', 18)
box off 
% saveas(gcf, fullfile(path, 'Suu_LES.png'))

return

figure
loglog(f*0.3/U{1}(11), S.*f/sqrt(uu{1}(11)), 'linewidth', 1.5, 'color', 'k')
ylim([1e-7 10])

% properties
ylabel('$$fS/\sqrt(\overline{u''^2})$$','Interpreter','latex'); 
xlabel('$$fB/U$$','Interpreter','latex');
set(gca, 'fontsize', 18)
box off
saveas(gcf, fullfile(path, 'Suu_LES.png'))

