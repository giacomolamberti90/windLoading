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

Tv_exp_mean = yLu_exp_mean./U_exp_mean;
Tv_exp_min  = yLu_exp_min./U_exp_max;
Tv_exp_max  = yLu_exp_max./U_exp_min;

zLu_exp_mean = dlmread(fullfile(path_EXP, 'zLuInlet'), '', [1  1  49  1]);
zLu_exp_min  = dlmread(fullfile(path_EXP, 'zLuInlet_min'), '', [1  1  49  1]);
zLu_exp_max  = dlmread(fullfile(path_EXP, 'zLuInlet_max'), '', [1  1  49  1]);

Tw_exp_mean = zLu_exp_mean./U_exp_mean;
Tw_exp_min  = zLu_exp_min./U_exp_max;
Tw_exp_max  = zLu_exp_max./U_exp_min;

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
path_SENS = '/home/storage/LES/highRise/fine/00deg/Smagorinsky/';
Nprobes = 21;
NLES = 27;
Nback = 5;
Nwindows = 8;

for num = 1:27%[5,11,13,14,15,17,23]
    
    sprintf('workdir = %i', num)
    path_LES = fullfile(path_SENS, sprintf('workdir.%i', num));
    
    U{num}  = dlmread(fullfile(path_LES, '/U_x_3p5m.out'));    
    
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
    
%     % mean velocity =======================================================
%     if num == 14
%         if floor(tprobes(i)) == tprobes(i)
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2d/UMean', tprobes(end))));
%         else
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2.2f/UMean', tprobes(end))));
%         end            
%     else
%         if floor(tprobes(i)) == tprobes(i)
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2d/UMean', tprobes(end))));
%         else
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2.1f/UMean', tprobes(end))));
%         end
%     end    
%     format  = repmat('%f', 1, (Nprobes+1)*3);
%     data_U_1p75 = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes+2);
%     fclose(fid);
%     
%     % reynolds stresses ===================================================
%     if num == 14
%         if floor(tprobes(i)) == tprobes(i)
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2d/UPrime2Mean', tprobes(end))));
%         else
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2.2f/UPrime2Mean', tprobes(end))));
%         end            
%     else
%         if floor(tprobes(i)) == tprobes(i)
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2d/UPrime2Mean', tprobes(end))));
%         else
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z1p75/%2.1f/UPrime2Mean', tprobes(end))));
%         end
%     end    
%     format  = repmat('%f', 1, (Nprobes+1)*6);
%     data_Ri_1p75 = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes+2);
%     fclose(fid);
%     
%     for i = 1:Nprobes 
%         % mean velocity
%         U_1p75{num}(i) = data_U_1p75{2+3*(i-1)}(end);
%         % velocity component
%         uu_1p75{num}(i) = data_Ri_1p75{2+6*(i-1)}(end);
%         vv_1p75{num}(i) = data_Ri_1p75{5+6*(i-1)}(end);
%         ww_1p75{num}(i) = data_Ri_1p75{7+6*(i-1)}(end);     
%         uv_1p75{num}(i) = data_Ri_1p75{3+6*(i-1)}(end);        
%     end
    
    % time scale ==========================================================       
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
        Tu_1p75{num}(i) = computeTimeScale(time-time(1), data_L_1p75(:,2+3*(i-1)), Nwindows);
        Tv_1p75{num}(i) = computeTimeScale(time-time(1), data_L_1p75(:,3+3*(i-1)), Nwindows);
        Tw_1p75{num}(i) = computeTimeScale(time-time(1), data_L_1p75(:,4+3*(i-1)), Nwindows);
    end
    clear data_Li_1p75 data_L_1p75 tprobes
    
    % ======================================================================================================================
    
    % time-scales =========================================================
    dire = dir(fullfile(path_LES, sprintf('/postProcessing/probes_z_1p45/')));
    for i = 1:length(dire)-2
        tprobes(i) = str2num(dire(2+i).name);
    end
    tprobes = sort(tprobes);

%     % mean velocity =======================================================
%     if num == 14
%         if floor(tprobes(i)) == tprobes(i)
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2d/UMean', tprobes(end))));
%         else
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2.2f/UMean', tprobes(end))));
%         end            
%     else
%         if floor(tprobes(i)) == tprobes(i)
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2d/UMean', tprobes(end))));
%         else
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2.1f/UMean', tprobes(end))));
%         end
%     end    
%     format  = repmat('%f', 1, (Nprobes+1)*3);
%     data_U_1p45 = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes+2);
%     fclose(fid);
%     
%     % reynolds stresses ===================================================
%     if num == 14
%         if floor(tprobes(i)) == tprobes(i)
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2d/UPrime2Mean', tprobes(end))));
%         else
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2.2f/UPrime2Mean', tprobes(end))));
%         end            
%     else
%         if floor(tprobes(i)) == tprobes(i)
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2d/UPrime2Mean', tprobes(end))));
%         else
%             fid = fopen(fullfile(path_LES, sprintf('postProcessing/probes_z_1p45/%2.1f/UPrime2Mean', tprobes(end))));
%         end
%     end    
%     format  = repmat('%f', 1, (Nprobes+1)*6);
%     data_Ri_1p45 = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes+2);
%     fclose(fid);
%     
%     for i = 1:Nprobes 
%         % mean velocity
%         U_1p45{num}(i) = data_U_1p45{2+3*(i-1)}(end);
%         % velocity component
%         uu_1p45{num}(i) = data_Ri_1p45{2+6*(i-1)}(end);
%         vv_1p45{num}(i) = data_Ri_1p45{5+6*(i-1)}(end);
%         ww_1p45{num}(i) = data_Ri_1p45{7+6*(i-1)}(end);     
%         uv_1p45{num}(i) = data_Ri_1p45{3+6*(i-1)}(end);        
%     end    
     
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
        Tu_1p45{num}(i) = computeTimeScale(time-time(1), data_L_1p45(:,2+3*(i-1)), Nwindows);
        Tv_1p45{num}(i) = computeTimeScale(time-time(1), data_L_1p45(:,3+3*(i-1)), Nwindows);
        Tw_1p45{num}(i) = computeTimeScale(time-time(1), data_L_1p45(:,4+3*(i-1)), Nwindows);
        %Lu_1p45{num}(i) = Tu_1p45{num}(i);
    end
    clear data_Li_1p45 data_L_1p45 tprobes

%     U{num}  = (U_1p45{num} + U_1p75{num})/2;    
%     
%     uu{num} = (uu_1p45{num} + uu_1p75{num})/2;
%     vv{num} = (vv_1p45{num} + vv_1p75{num})/2;
%     ww{num} = (ww_1p45{num} + ww_1p75{num})/2;
%     
%     Iu{num} = sqrt(uu{num})./U{num};
%     Iv{num} = sqrt(vv{num})./U{num};
%     Iw{num} = sqrt(ww{num})./U{num};
%     
%     TKE{num} = 0.5 * (uu{num} + vv{num} + ww{num});    
    
    Tu{num} = (Tu_1p45{num}' + Tu_1p75{num}')/2;
    Tv{num} = (Tv_1p45{num}' + Tv_1p75{num}')/2;
    Tw{num} = (Tw_1p45{num}' + Tw_1p75{num}')/2;    
end

%%

effect = 'z_0';

for p = 4:6
            
    figure
    hold on
    set(gca, 'fontsize', 22)

    if p == 1
        variable = uu;
        var_exp_mean = uu_exp_mean;
        var_exp_min = uu_exp_min;
        var_exp_max = uu_exp_max;
        xlabel('$$\overline{u''^2} [m/s]$$', 'interpreter', 'latex');
    elseif p == 2
        variable = vv;
        var_exp_mean = vv_exp_mean;
        var_exp_min = vv_exp_min;
        var_exp_max = vv_exp_max;
        xlabel('$$\overline{v''^2} [m/s]$$', 'interpreter', 'latex');
    elseif p == 3
        variable = ww;
        var_exp_mean = ww_exp_mean;
        var_exp_min = ww_exp_min;
        var_exp_max = ww_exp_max;
        xlabel('$$\overline{w''^2} [m/s]$$', 'interpreter', 'latex');
    elseif p == 4
        variable = Tu;
        var_exp_mean = Tu_exp_mean;
        var_exp_min = Tu_exp_min;
        var_exp_max = Tu_exp_max;
        xlabel('$$T_u [s]$$', 'interpreter', 'latex');
    elseif p == 5
        variable = Tv;
        var_exp_mean = Tv_exp_mean;
        var_exp_min = Tv_exp_min;
        var_exp_max = Tv_exp_max;
        xlabel('$$T_v [s]$$', 'interpreter', 'latex');
    elseif p == 6
        variable = Tw;
        var_exp_mean = Tw_exp_mean;
        var_exp_min = Tw_exp_min;
        var_exp_max = Tw_exp_max;
        xlabel('$$T_w [s]$$', 'interpreter', 'latex');
    end 
    
    var_min  = [];
    var_mean = [];
    var_max  = [];
    
    if effect == 'z_0'    
        for i = 1:9
            index = [3*(i-1) + 1, 3*(i-1) + 2, 3*(i-1) + 3];
            var_min  = [var_min;  variable{index(1)}'];
            var_mean = [var_mean; variable{index(2)}'];
            var_max  = [var_max;  variable{index(3)}'];
        end           
    elseif effect == 'tke'
        for i = 1:3
            for j = 1:3
                index = [i + (j-1)*9, i+3 + (j-1)*9, i+6 + (j-1)*9];
                var_min  = [var_min;  variable{index(1)}'];
                var_mean = [var_mean; variable{index(2)}'];
                var_max  = [var_max;  variable{index(3)}'];
            end
        end                
    elseif effect == 'T_u'
        for i = 1:9
            index = [i, i+9, i+18];
            var_min  = [var_min;  variable{index(1)}'];
            var_mean = [var_mean; variable{index(2)}'];
            var_max  = [var_max;  variable{index(3)}'];
        end
    end
    
    % end
    plot(mean(var_min),  linspace(0,4,21), 'b', 'linewidth', 2)
    plot(mean(var_mean), linspace(0,4,21), 'k', 'linewidth', 2)
    plot(mean(var_max),  linspace(0,4,21), 'r', 'linewidth', 2)

    plot(var_exp_mean, z_exp, '--k')
    y = [z_exp; flipud(z_exp)];
    x = [var_exp_min; flipud(var_exp_max)];
    fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
    
    if effect == 'z_0'
        legend({'$$z_0$$ - min','$$z_0$$ - mean','$$z_0$$ - max'}, 'interpreter', 'latex')
    elseif effect == 'tke'
        legend({'$$k$$ - min','$$k$$ - mean','$$k$$ - max'}, 'interpreter', 'latex')
    elseif effect == 'T_u'
        legend({'$$T_u$$ - min','$$T_u$$ - mean','$$T_u$$ - max'}, 'interpreter', 'latex')
    end
    
    legend boxoff
    ylabel('$$z [m]$$', 'interpreter', 'latex')    
    ylim([0 2])
end

return

%% figure -----------------------------------------------------------------
path = '/home/giacomol/Desktop/Research/slides/10-15-19/figures';

z = linspace(0,4);

figure
for i = 1:9
    index = [3*(i-1) + 1, 3*(i-1) + 2, 3*(i-1) + 3];

    subplot(3,3,i)    
    hold on
    set(gca, 'fontsize', 18)
    for num = LES
        plot(U{num}, z, 'linewidth', 2)
    end
    plot(U_exp_mean, z_exp, '--k')
    y = [z_exp; flipud(z_exp)];
    x = [U_exp_min; flipud(U_exp_max)];
    fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
    if i == 1
        legend('min','mean','max')
    end
    ylabel('$$z [m]$$', 'interpreter', 'latex')
    xlabel('$$\overline{u} [m/s]$$', 'interpreter', 'latex');
    ylim([0 2])
    % saveas(gcf, fullfile(path, 'U_k.png'))
end


figure
k = 1;
for i = 1:3
    for j = 1:3
        LES = [i + (j-1)*9, i+3 + (j-1)*9, i+6 + (j-1)*9];

        subplot(3,3,k)
        hold on
        set(gca, 'fontsize', 18)
        for num = LES
            plot(TKE{num}, z, 'linewidth', 2)
        end
        plot(k_exp_mean, z_exp, '--k')
        y = [z_exp; flipud(z_exp)];
        x = [k_exp_min; flipud(k_exp_max)];
        fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
        if i == 1
            legend('min','mean','max')
        end    
        ylabel('$$z [m]$$', 'interpreter', 'latex')
        xlabel('$$k [m^2/s^2]$$', 'interpreter', 'latex');
        axis([0 2 0 2])
        % saveas(gcf, fullfile(path, 'k_k.png'))
        k = k+1;
    end
end

% Tu
z = linspace(0,4,21);

figure
for i = 1:9
    LES = [i, i+9, i+18];

    subplot(3,3,i)    
    hold on
    set(gca, 'fontsize', 18)
    for num = LES
        plot(Tu{num}, z, 'linewidth', 2)
    end
    plot(Tu_exp_mean, z_exp, '--k')
    y = [z_exp; flipud(z_exp)];
    x = [Tu_exp_min; flipud(Tu_exp_max)];
    fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
    if i == 1
        legend('min','mean','max')
    end    
    ylabel('$$z [m]$$', 'interpreter', 'latex')
    xlabel('$$T_u [s]$$', 'interpreter', 'latex');
    ylim([0 2])
    % saveas(gcf, fullfile(path, 'Tu_k.png'))
end

return

%%
z = linspace(0,4);

figure
hold on
set(gca, 'fontsize', 18)
for num = 1:27
    plot(U{num}, z, 'r')
end
plot(U{14}, z, '-k', 'linewidth', 2)
plot(U_exp_mean, z_exp, '--k')
y = [z_exp; flipud(z_exp)];
x = [U_exp_min; flipud(U_exp_max)];
fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$U [m/s]$$', 'interpreter', 'latex');
ylim([0 2])

figure
hold on
set(gca, 'fontsize', 18)
for num = 1:27
    plot(TKE{num}, z, 'r')
end
plot(TKE{14}, z, '-k', 'linewidth', 2)
plot(k_exp_mean, z_exp, '--k')
y = [z_exp; flipud(z_exp)];
x = [k_exp_min; flipud(k_exp_max)];
fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$k [m^2/s^2]$$', 'interpreter', 'latex');
ylim([0 2])

z = linspace(0,4,21);

figure
hold on
set(gca, 'fontsize', 18)
for num = 1:27
    plot(Tu{num}, z, 'r')
end
plot(Tu{14}, z, '-k', 'linewidth', 2)
plot(Tu_exp_mean, z_exp, '--k')
y = [z_exp; flipud(z_exp)];
x = [Tu_exp_min; flipud(Tu_exp_max)];
fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$T_u [s]$$', 'interpreter', 'latex');
ylim([0 2])

return

% R11
figure
subplot(1,3,1)
hold on
set(gca, 'fontsize', 18)
for num = 1:27
    plot(Iu{num}, z, 'r')
end
plot(Iu_exp_mean, z_exp, '--k')
y = [z_exp; flipud(z_exp)];
x = [Iu_exp_min; flipud(Iu_exp_max)];
fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$I_u$$', 'interpreter', 'latex');
ylim([0 2])

% R22
subplot(1,3,2)
hold on
set(gca, 'fontsize', 18)
for num = 1:27
    plot(Iv{num}, z, 'r')
end
plot(Iv_exp_mean, z_exp, '--k')
y = [z_exp; flipud(z_exp)];
x = [Iv_exp_min; flipud(Iv_exp_max)];
fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$I_v$$', 'interpreter', 'latex');
ylim([0 2])

% R33
subplot(1,3,3)
hold on
set(gca, 'fontsize', 18)
for num = 1:27
    plot(Iw{num}, z, 'r')
end
plot(Iw_exp_mean, z_exp, '--k')
y = [z_exp; flipud(z_exp)];
x = [Iw_exp_min; flipud(Iw_exp_max)];
fill(x, y, 'k', 'linestyle', 'none', 'facealpha', .4);
ylabel('$$z [m]$$', 'interpreter', 'latex')
xlabel('$$I_w$$', 'interpreter', 'latex');
ylim([0 2])
