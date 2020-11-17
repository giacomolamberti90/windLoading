function convergenceLES(path_LES, ang, N, nT)

    % check convergence every t_* = 10 (~2sec)
    %                         t_* = Tu_@2m*U_@2m/H = 0.4
            
%% data
% read probes
dt = 0.0001;
Nprobes_V = 7; NV = 1;

data_A = [];
data_B = [];

if ang == 0 || ang == 180
    Nprobes_A = 224;
    Nprobes_B = 223;
elseif ang == 20
    Nprobes_A = 1433;
elseif ang == 50;
    Nprobes_A = 2023;
elseif ang == 60;
    Nprobes_A = 2132;
elseif ang == 70;
    Nprobes_A = 2253;    
elseif ang == 80;
    Nprobes_A = 2368;
elseif ang == 90;
    Nprobes_A = 2394;    
end

% time in probes
dire = dir(fullfile(path_LES, sprintf('/postProcessing/A%d', ang)));
for i = 1:length(dire)-2
    tprobes(i) = str2num(dire(2+i).name);
end
tprobes = sort(tprobes);

% reference velocity
fid    = fopen(fullfile(path_LES, sprintf('postProcessing/probes_x_2m/%s/UMean', num2str(tprobes(end)))));
format = repmat('%f', 1, 3*(Nprobes_V+1));
data_U = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_V+2);
fclose(fid); 

U = data_U{2+3*(NV-1)}(end);
sprintf('Uref = %f', U)

%% 0 deg - mean
% read mean and variance and time series
for i = 1:length(tprobes)
    % time series
    fid    = fopen(fullfile(path_LES, sprintf('postProcessing/A%d/%s/pMean', ang, num2str(tprobes(i)))));
    format = repmat('%f', 1, Nprobes_A+1);
    data   = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_A+2);
    fclose(fid); 

    data_A = [data_A; cell2mat(data)];
    clear data
    
    if ang == 0 || ang == 180
        fid    = fopen(fullfile(path_LES, sprintf('postProcessing/B%d/%s/pMean', ang, num2str(tprobes(i)))));
        format = repmat('%f', 1, Nprobes_B+1);
        data   = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_B+2);
        fclose(fid);

        data_B = [data_B; cell2mat(data)];
        clear data
    end
end  

if ang == 0 || ang == 180
    cp_A = data_A(:,2:end)/(0.5*U^2); %cumsum(data_A(:,2:end))./repmat([1:size(data_A,1)]',1,224)/(0.5*U^2);
    cp_B = data_B(:,2:end)/(0.5*U^2); %cumsum(data_B(:,2:end))./repmat([1:size(data_B,1)]',1,223)/(0.5*U^2);

    [~, ~, ind_A] = intersect('A1313', taps_A_ref, 'rows', 'stable');
    [~, ~, ind_B] = intersect('B0808', taps_B_ref, 'rows', 'stable');

    err_A = max(abs(cp_A(end,ind_A)-cp_A(end-nT/dt,ind_A))./abs(cp_A(end,ind_A)));
    err_B = max(abs(cp_B(end,ind_B)-cp_B(end-nT/dt,ind_B))./abs(cp_B(end,ind_B)));
else
    cp_A = data_A(:,2:end)/(0.5*U^2);
%     ind_A = 1030;
%     err_A = max(abs(cp_A(end,ind_A)-cp_A(end-nT/dt,ind_A))./abs(cp_A(end,ind_A)));
%     err_B = err_A;
end
% sprintf('The maximum difference in %6.4fs in the mean is %6.4f%%', data_A(end,1)-data_A(end-nT/dt,1), 100*max(err_A, err_B))

figure(1)
hold on
xlabel('$$t [s]$$', 'interpreter', 'latex')
ylabel('$$C_{p,mean}$$', 'interpreter', 'latex')
if ang == 0 || ang == 180
    plot(data_A(:,1), cp_A(:,ind_A))
    plot(data_B(:,1), cp_B(:,ind_B))
else
    plot(data_A(:,1), cp_A)
end
set(gca, 'fontsize', 18)
ylim([-1.5 0])

clear data_A data_B

%% 0 deg - rms
data_A = [];
data_B = [];

% read mean and variance and time series
for i = 1:length(tprobes)
    % time series
    fid    = fopen(fullfile(path_LES, sprintf('postProcessing/A%d/%s/pPrime2Mean', ang, num2str(tprobes(i)))));
    format = repmat('%f', 1, Nprobes_A+1);
    data   = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_A+2);
    fclose(fid); 

    data_A = [data_A; cell2mat(data)];
    clear data
    
    if ang == 0 || ang == 180
        fid    = fopen(fullfile(path_LES, sprintf('postProcessing/B%d/%s/pPrime2Mean', ang, num2str(tprobes(i)))));
        format = repmat('%f', 1, Nprobes_B+1);
        data   = textscan(fid, format, 'Delimiter','() ','MultipleDelimsAsOne', true, 'Headerlines', Nprobes_B+2);
        fclose(fid);

        data_B = [data_B; cell2mat(data)];
        clear data
    end
end  

if ang == 0 || ang == 180
    cp_A = data_A(:,2:end)/(0.5*U^2); %cumsum(data_A(:,2:end))./repmat([1:size(data_A,1)]',1,224)/(0.5*U^2);
    cp_B = data_B(:,2:end)/(0.5*U^2); %cumsum(data_B(:,2:end))./repmat([1:size(data_B,1)]',1,223)/(0.5*U^2);

    [~, ~, ind_A] = intersect('A1313', taps_A_ref, 'rows', 'stable');
    [~, ~, ind_B] = intersect('B0808', taps_B_ref, 'rows', 'stable');

    err_A = max(abs(cp_A(end,ind_A)-cp_A(end-nT/dt,ind_A))./abs(cp_A(end,ind_A)));
    err_B = max(abs(cp_B(end,ind_B)-cp_B(end-nT/dt,ind_B))./abs(cp_B(end,ind_B)));
else
    cp_A = data_A(:,2:end)/(0.5*U^2);
%     ind_A = 1030;
%     err_A = max(abs(cp_A(end,ind_A)-cp_A(end-nT/dt,ind_A))./abs(cp_A(end,ind_A)));
%     err_B = err_A;
end

% sprintf('The maximum difference in %6.4fs in the rms is %6.4f%%', data_A(end,1)-data_A(end-nT/dt,1), 100*max(err_A, err_B))

figure(2)
hold on
xlabel('$$t [s]$$', 'interpreter', 'latex')
ylabel('$$C_{p,rms}$$', 'interpreter', 'latex')
if ang == 0 || ang == 180
    plot(data_A(:,1), cp_A(:,ind_A))
    plot(data_B(:,1), cp_B(:,ind_B))
else
    plot(data_A(:,1), cp_A)
end
set(gca, 'fontsize', 18)
ylim([0 0.75])
