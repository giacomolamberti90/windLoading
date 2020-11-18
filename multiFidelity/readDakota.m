%% Read Dakota output

clear; close all; clc

Ntaps = 8384;
Nlevels = 4001;
Norder = 12;

filename = '../7order/myLES/skip5/run_mean_RANS.out';
run = 'rans7';

%% mean and std
[fl, p] = grep('-u -n', {'expansion:'}, filename);
lines = p.result;
for i = 2:size(lines,1)-1
   
    [firstcolumns, pos] = textscan(lines(i,:),'%s', 2); %first four columns
    cp = textscan(lines(i,pos+1:end), '%s'); %remaining numbers
    cp_mean(i-1) = str2num(char(cp{1}(1)));
    cp_std(i-1) = str2num(char(cp{1}(2)));
    clear data
end

fid = fopen(filename, 'rt');
s = textscan(fid, '%s', 'delimiter', '\n');

A = [cp_mean', cp_mean' - 1.96 * cp_std', cp_mean' + 1.96 * cp_std'];

filename = '../7order/myLES/skip5/PCE_mean_RANS_PCE7_stdCI.out';
dlmwrite(filename, A, 'delimiter', ' ');

%% CDF
tic
if run == 'multi'
    % find row id
    idx = find(strcmp(s{1}, 'Multifidelity UQ: approximated high fidelity results'), 1, 'first');
    idx1 = find(strcmp(s{1}(idx:end), 'Cumulative Distribution Function (CDF) for response_fn_1:'), 1, 'first');
    idx2 = find(strcmp(s{1}(idx:end), ['Cumulative Distribution Function (CDF) for response_fn_', num2str(Ntaps), ':']), 1, 'first');

    % gather text between indices
    AA = s{1}(idx + idx1 : idx + idx2 + Nlevels + 3);
    
    for i = 1:Ntaps
        i
        for j = 1:Nlevels-1
            x(j,i) = str2num(AA{(i-1)*Nlevels + 3*i + j}(1:18));
            f(j,i) = str2num(AA{(i-1)*Nlevels + 3*i + j}(19:end));
        end
        [f_ind, ind]  = unique(f(:,i));
        CI(i,:) = interp1(f_ind, x(ind, i), [0.025 0.975]);
        clear f_ind ind
    end
    
else
    % find row id
    idx1 = find(strcmp(s{1}, 'Cumulative Distribution Function (CDF) for response_fn_1:'), 1, 'first');
    idx2 = find(strcmp(s{1}, ['Cumulative Distribution Function (CDF) for response_fn_', num2str(Ntaps), ':']), 1, 'first');

    % gather text between indices
    AA = s{1}(idx1 : idx2 + Nlevels + 3);  
    
    for i = 1:Ntaps
        i
        for j = 1:Nlevels
            x(j,i) = str2num(AA{(i-1)*Nlevels + 3*i + j}(1:18));
            f(j,i) = str2num(AA{(i-1)*Nlevels + 3*i + j}(19:end));
        end
        [f_ind, ind]  = unique(f(:,i));
        CI(i,:) = interp1(f_ind, x(ind, i), [0.025 0.975]);
        clear f_ind ind
    end
end
toc

A = [cp_mean', CI];

filename = '../7order/myLES/skip5/PCE_mean_RANS_PCE7.out';
dlmwrite(filename, A, 'delimiter', ' ');

%% Build PCE
PCE = @(order, angle) legendreP(order, angle);

angles = [0, 10, 20, 40, 45, 60, 80, 90];

% angles = [0.2775616554, 1.777892908, 5.019334521, 10.14314988, 17.0503674,...
%           25.45903128, 34.94759911, 45, 55.05240089, 64.54096872, 72.9496326,...
%           79.85685012, 84.98066548, 88.22210709, 89.72243834];

left = 0;
right = 90;

idx3 = find(strcmp(s{1}, 'Coefficients of Polynomial Chaos Expansion for response_fn_1:'), 1, 'first');
idx4 = find(strcmp(s{1}, ['Coefficients of Polynomial Chaos Expansion for response_fn_',num2str(Ntaps),':']),1,'first');

BB = s{1}(idx3 : idx4 + Norder + 3);

for i = 1:Ntaps
    for j = 1:Norder
        f_B(j,i) = str2num(BB{(i-1)*Norder + 3*i + j}(1:18));
    end
end

Cp_mean_angle = zeros(length(angles), Ntaps);

for i = 1:length(angles)
    tic
    angles(i)
    parfor j = 1:Ntaps
        val = 0;
        for k = 1:Norder;
            ang = (2*angles(i)-left-right) / (right-left);
            val = val + (PCE(k-1, ang) * f_B(k,j));
        end
        Cp_mean_angle(i,j) = val;
    end
    toc
end

dlmwrite('../7order/myLES/skip5/PCE_mean_RANS_PCE7_angles.txt', Cp_mean_angle', 'delimiter', ' ')
