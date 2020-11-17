function [tileA, tileB] = pressureUQ(angle, stat)

    % Function that read UQ files RANS of tiles A-B and returns:
        
        %structure containing:
        %        coords: coordinates of pressure taps           [Ntaps,2]
        %          mean: meanUQ of pressure coefficients        [Ntaps,1]
        %           var: varianceUQ of pressure coefficients    [Ntaps,1]
        
%% taps list - reference
tapsList   = load('/home/giacomol/Desktop/Research/windLoading/windTunnel/matlabCode/ArupGUI3.0/Orifizi-v2.mat');
taps_A_ref = cell2mat(tapsList.nomi{1}');
taps_B_ref = cell2mat(tapsList.nomi{2}');
taps_D_ref = cell2mat(tapsList.nomi{4}');

%% data
Ntaps = 224;
path  = '/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/';
if angle == 0
    pathUQ = '/home/giacomol/Desktop/Research/windLoading/RANS/highRise/UQ/00deg/aleatoric';    
    % taps coordinates
    coords_A = load(fullfile(path, 'coords_A0'));
    coords_B = load(fullfile(path, 'coords_B0'));
    % files
    if stat == 'mean'
        fileA = 'run_A0.out';
        fileB = 'run_B0.out';
        Nlevels = 751;
    elseif stat == 'rms1'
        fileA = 'run_rms_PH_A0.out';
        fileB = 'run_rms_PH_B0.out';   
        Nlevels = 501;
    elseif stat == 'rms2'
        fileA = 'run_rms_S_A0.out';
        fileB = 'run_rms_S_B0.out';   
        Nlevels = 501;        
    elseif stat == 'rms3'
        fileA = 'run_rms_P_A0.out';
        fileB = 'run_rms_P_B0.out';   
        Nlevels = 501;
    elseif stat == 'rms4'
        fileA = 'run_rms_RW_A0.out';
        fileB = 'run_rms_RW_B0.out';   
        Nlevels = 501;     
    end
elseif angle == 180  
    pathUQ = '/home/giacomol/Desktop/Research/windLoading/RANS/highRise/UQ/00deg/aleatoric';    
    % taps coordinates
    coords_A = load(fullfile(path, 'coords_A180'));
    coords_B = load(fullfile(path, 'coords_B180'));
    % files
    if stat == 'mean'
        fileA = 'run_A180.out';
        fileB = 'run_B180.out';
        Nlevels = 801;
    elseif stat == 'rms1'
        fileA = 'run_rms_PH_A180.out';
        fileB = 'run_rms_PH_B180.out';   
        Nlevels = 501;
    elseif stat == 'rms2'
        fileA = 'run_rms_S_A180.out';
        fileB = 'run_rms_S_B180.out';   
        Nlevels = 501;        
    elseif stat == 'rms3'
        fileA = 'run_rms_P_A180.out';
        fileB = 'run_rms_P_B180.out';   
        Nlevels = 501;
    elseif stat == 'rms4'
        fileA = 'run_rms_RW_A180.out';
        fileB = 'run_rms_RW_B180.out';   
        Nlevels = 501;     
    end    
elseif angle == 20
    pathUQ = '/home/giacomol/Desktop/Research/windLoading/RANS/highRise/UQ/20deg/aleatoric';    
    % taps coordinates
    coords_A = load(fullfile(path, 'coords_A0'));
    coords_B = load(fullfile(path, 'coords_B0'));
    % files
    fileA = 'run_A20.out';
    fileB = 'run_B20.out';
    Nlevels = 800;    
else
    disp('RANS only at 0,20,180 [deg]')
    tileA = [];
    tileB = [];
    return
end
pdata_A = fullfile(pathUQ, fileA);
pdata_B = fullfile(pathUQ, fileB);

%% UQ
% mean and std
[fl, p] = grep('-u -n', {'expansion:'}, pdata_A);
lines = p.result;
for i = 2:size(lines,1)-1
   
    [firstcolumns, pos] = textscan(lines(i,:),'%s', 2);    %first four columns
    p_A           = textscan(lines(i,pos+1:end), '%s');   %remaining numbers
    cpmean_A(i-1) = str2num(char(p_A{1}(1)));
    clear data
end

% mean and std
[fl, p] = grep('-u -n', {'expansion:'}, pdata_B);
lines = p.result;
for i = 2:size(lines,1)-1
   
    [firstcolumns, pos] = textscan(lines(i,:),'%s', 2); %first four columns
    p_B           = textscan(lines(i,pos+1:end), '%s');   %remaining numbers
    cpmean_B(i-1) = str2num(char(p_B{1}(1)));
    clear data
end

%% CDF

% UQ output file - tile A
fid = fopen(pdata_A, 'rt');
s   = textscan(fid, '%s', 'delimiter', '\n');

% find row id
idx1 = find(strcmp(s{1}, 'Cumulative Distribution Function (CDF) for response_fn_1:'), 1, 'first');
idx2 = find(strcmp(s{1}, ['Cumulative Distribution Function (CDF) for response_fn_',num2str(Ntaps),':']),1,'first');

% gather text between indices
AA = s{1}(idx1:idx2+Nlevels+3);

for i = 1:Ntaps
    for j = 1:Nlevels
        x_A(j,i) = str2num(AA{(i-1)*Nlevels+3*i+j}(1:18));
        f_A(j,i) = str2num(AA{(i-1)*Nlevels+3*i+j}(19:end));
    end
    [f, ind]  = unique(f_A(:,i));
    CI_A(i,:) = interp1(f, x_A(ind,i), [0.025 0.975]);
    clear f ind
end

% UQ output file - tile B
fid = fopen(pdata_B, 'rt');
s   = textscan(fid, '%s', 'delimiter', '\n');

% find row id
idx1 = find(strcmp(s{1}, 'Cumulative Distribution Function (CDF) for response_fn_1:'), 1, 'first');
idx2 = find(strcmp(s{1}, ['Cumulative Distribution Function (CDF) for response_fn_',num2str(Ntaps-1),':']),1,'first');

% gather text between indices
AA = s{1}(idx1:idx2+Nlevels+3);

for i = 1:Ntaps-1
    for j = 1:Nlevels
        x_B(j,i) = str2num(AA{(i-1)*Nlevels+3*i+j}(1:18));
        f_B(j,i) = str2num(AA{(i-1)*Nlevels+3*i+j}(19:end));
    end
    [f, ind]  = unique(f_B(:,i));
    CI_B(i,:) = interp1(f, x_B(ind,i), [0.025 0.975]);
    clear f ind
end

%% Sobol indices
% UQ output file - tile A
fid = fopen(pdata_A, 'rt');
s   = textscan(fid, '%s', 'delimiter', '\n');

% find row id
idx1 = find(strcmp(s{1}, 'response_fn_1 Sobol'' indices:'), 1, 'first');
idx2 = find(strcmp(s{1}, ['response_fn_',num2str(Ntaps), ' Sobol'' indices:']), 1, 'first');

% gather text between indices
AA = s{1}(idx1:idx2+9);

k  = 1;
ii = 1;
for i = 1:Ntaps
    
    % check if Sobol indices are computed
    xx = strmatch(['response_fn_', num2str(k),' Sobol'' indices not available due to negligible variance'], AA{ii});
    if xx == 1
        
        FN.sobol.x1(k)     = 0;
        FN.sobol.x2(k)     = 0;
        FN.sobol.x2(k)     = 0;
        FN.sobol.x1x2(k)   = 0;
        FN.sobol.x1x3(k)   = 0;
        FN.sobol.x2x3(k)   = 0;
        FN.sobol.x1x2x3(k) = 0;
        
        ii = ii+1;
    else
        
        AA_num = str2num(AA{ii+2}(1:end-3));
        if AA_num > 1
            xxx
        end
    
        FN.sobol.x1(k)=AA_num(1);
        AA_num = str2num(AA{ii+3}(1:end-3));
        if AA_num > 1
            xxx
        end
    
        FN.sobol.x2(k) = AA_num(1);
    
        AA_num         = str2num(AA{ii+4}(1:end-4));
        FN.sobol.x3(k) = AA_num(1);
     
        AA_num           = str2num(AA{ii+6}(1:end-6));
        FN.sobol.x1x2(k) = AA_num(1);
    
        AA_num           = str2num(AA{ii+7}(1:end-6));
        FN.sobol.x1x3(k) = AA_num(1);
    
        AA_num           = str2num(AA{ii+8}(1:end-6));
        FN.sobol.x2x3(k) = AA_num(1);
    
        AA_num             = str2num(AA{ii+9}(1:end-9));
        FN.sobol.x1x2x3(k) = AA_num(1);
    
    
        ii = ii+10;
    end
    k = k+1;
end

Si_A(:,1) = FN.sobol.x1;
Si_A(:,2) = FN.sobol.x2;
Si_A(:,3) = FN.sobol.x3;
Si_A(:,4) = FN.sobol.x1x2;
Si_A(:,5) = FN.sobol.x2x3;
Si_A(:,6) = FN.sobol.x1x3;

clear FN

% UQ output file - tile A
fid = fopen(pdata_B, 'rt');
s   = textscan(fid, '%s', 'delimiter', '\n');

% find row id
idx1 = find(strcmp(s{1}, 'response_fn_1 Sobol'' indices:'), 1, 'first');
idx2 = find(strcmp(s{1}, ['response_fn_',num2str(Ntaps-1), ' Sobol'' indices:']), 1, 'first');

% gather text between indices
AA = s{1}(idx1:idx2+9);

k  = 1;
ii = 1;
for i = 1:Ntaps-1
    
    % check if Sobol indices are computed
    xx = strmatch(['response_fn_', num2str(k),' Sobol'' indices not available due to negligible variance'], AA{ii});
    if xx == 1
        
        FN.sobol.x1(k)     = 0;
        FN.sobol.x2(k)     = 0;
        FN.sobol.x2(k)     = 0;
        FN.sobol.x1x2(k)   = 0;
        FN.sobol.x1x3(k)   = 0;
        FN.sobol.x2x3(k)   = 0;
        FN.sobol.x1x2x3(k) = 0;
        
        ii = ii+1;
    else
        
        AA_num = str2num(AA{ii+2}(1:end-3));
        if AA_num > 1
            xxx
        end
    
        FN.sobol.x1(k)=AA_num(1);
        AA_num = str2num(AA{ii+3}(1:end-3));
        if AA_num > 1
            xxx
        end
    
        FN.sobol.x2(k) = AA_num(1);
    
        AA_num         = str2num(AA{ii+4}(1:end-4));
        FN.sobol.x3(k) = AA_num(1);
     
        AA_num           = str2num(AA{ii+6}(1:end-6));
        FN.sobol.x1x2(k) = AA_num(1);
    
        AA_num           = str2num(AA{ii+7}(1:end-6));
        FN.sobol.x1x3(k) = AA_num(1);
    
        AA_num           = str2num(AA{ii+8}(1:end-6));
        FN.sobol.x2x3(k) = AA_num(1);
    
        AA_num             = str2num(AA{ii+9}(1:end-9));
        FN.sobol.x1x2x3(k) = AA_num(1);
    
    
        ii = ii+10;
    end
    k = k+1;
end

Si_B(:,1) = FN.sobol.x1;
Si_B(:,2) = FN.sobol.x2;
Si_B(:,3) = FN.sobol.x3;
Si_B(:,4) = FN.sobol.x1x2;
Si_B(:,5) = FN.sobol.x2x3;
Si_B(:,6) = FN.sobol.x1x3;

%% define structures
% tileA
tileA.CI        = 'on';
tileA.CI95_mean = CI_A;
tileA.taps      = taps_A_ref;
tileA.coords    = coords_A;
tileA.mean      = cpmean_A';
tileA.sobol     = Si_A;
tileA.xcdf      = x_A;
tileA.fcdf      = f_A;
%tileB
tileB.CI        = 'on';
tileB.CI95_mean = CI_B;
tileB.taps      = taps_B_ref;
tileB.coords    = coords_B;
tileB.mean      = cpmean_B';
tileB.sobol     = Si_B;
tileB.xcdf      = x_B;
tileB.fcdf      = f_B;
 