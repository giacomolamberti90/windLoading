function highRise = pressureUQ_all(angle, comp)

    % Function that given angle and component returns:
        
        %structure containing:
        %        coords: coordinates of pressure taps           [Ntaps,2]
        %          mean: meanUQ of pressure coefficients        [Ntaps,1]
        %           var: varianceUQ of pressure coefficients    [Ntaps,1]
        
%% taps list - reference
tapsList   = load('/home/giacomol/Desktop/Research/windLoading/windTunnel/matlabCode/ArupGUI3.0/Orifizi-v2.mat');
taps_A_ref = cell2mat(tapsList.nomi{1}');
taps_B_ref = cell2mat(tapsList.nomi{2}');
taps_C_ref = cell2mat(tapsList.nomi{3}');

%% data
path  = '/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/';
if angle == 0
    pathUQ = '/home/giacomol/Desktop/Research/windLoading/RANS/highRise/UQ/00deg/';
    
    % taps coordinates
    coords = load(fullfile(path, 'coords_all'));
    
    % files
    file = 'run_all.out';
    if comp == 'bas'
       Nlevels = 2001;
    elseif comp == 'c1c'
       file = 'run_all_delta05.out'; 
       Nlevels = 1001;
    else
       Nlevels = 1001;
    end
    
elseif angle == 20
    pathUQ = '/home/giacomol/Desktop/Research/windLoading/RANS/highRise/UQ/20deg/';
    % taps coordinates
    coords = load(fullfile(path, 'coords_all'));
    
    % files
    file = 'run_all.out';
    if comp == 'c1c'
        file = 'run_all_delta05.out';
        Nlevels = 1001;
    else
        Nlevels = 1001;
    end
    
else
    disp('RANS only at 0,20 [deg]')
    highRise = [];
    return
end

if comp == 'bas'
    pathUQ = fullfile(pathUQ, sprintf('aleatoric/'));
else
    pathUQ = fullfile(pathUQ, sprintf('epistemic/%s/', comp));
    pathUQ = fullfile(pathUQ, 'delta0p2/');
end

pdata = fullfile(pathUQ, file);

%% UQ
% mean and std
[fl, p] = grep('-u -n', {'expansion:'}, pdata);
lines = p.result;
for i = 2:size(lines,1)-1
   
    [firstcolumns, pos] = textscan(lines(i,:),'%s', 2);       %first four columns
    pressure            = textscan(lines(i,pos+1:end), '%s'); %remaining numbers
    cpmean(i-1)         = str2num(char(pressure{1}(1)));
    cpstd(i-1)          = str2num(char(pressure{1}(2)));
    clear data
end

%% CDF
% UQ output file - tile A
fid = fopen(pdata, 'rt');
s   = textscan(fid, '%s', 'delimiter', '\n');

% find row id
idx1 = find(strcmp(s{1}, 'Cumulative Distribution Function (CDF) for response_fn_1:'), 1, 'first');
idx2 = find(strcmp(s{1}, 'Cumulative Distribution Function (CDF) for response_fn_1530:'),1,'first');

% gather text between indices
AA = s{1}(idx1:idx2+Nlevels+3);

for i = 1:size(coords,1)
    for j = 1:Nlevels
        x_pdf(j,i) = str2num(AA{(i-1)*Nlevels+3*i+j}(1:18));
        f_pdf(j,i) = str2num(AA{(i-1)*Nlevels+3*i+j}(19:end));
    end
    [funique, ind]  = unique(f_pdf(:,i));
    CI(i,:) = interp1(funique, x_pdf(ind,i), [0.025 0.975]);
    clear funique ind
end

%% Sobol indices
% UQ output file - tile A
fid = fopen(pdata, 'rt');
s   = textscan(fid, '%s', 'delimiter', '\n');

% find row id
idx1 = find(strcmp(s{1}, 'response_fn_1 Sobol'' indices:'), 1, 'first');
idx2 = find(strcmp(s{1}, 'response_fn_1530 Sobol'' indices:'), 1, 'first');

% gather text between indices
AA = s{1}(idx1:idx2+9);

k  = 1;
ii = 1;
for i = 1:size(coords,1)
    
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

Sobol(:,1) = FN.sobol.x1;
Sobol(:,2) = FN.sobol.x2;
Sobol(:,3) = FN.sobol.x3;
Sobol(:,4) = FN.sobol.x1x2;
Sobol(:,5) = FN.sobol.x2x3;
Sobol(:,6) = FN.sobol.x1x3;

%% define structures
N = size(coords,1) - 2*(length(taps_A_ref) + length(taps_B_ref) + length(taps_C_ref));

taps_all = [taps_A_ref; taps_A_ref; taps_B_ref; taps_B_ref; taps_C_ref; taps_C_ref];
taps_all = [taps_all; repmat('XXXXX', N, 1)];

% tileA
highRise.CI        = 'on';
highRise.CI95_mean = CI;
highRise.taps      = taps_all;
highRise.coords    = coords;
highRise.mean      = cpmean';
highRise.std       = cpstd';
highRise.sobol     = Sobol;
highRise.xcdf      = x_pdf;
highRise.fcdf      = f_pdf;
