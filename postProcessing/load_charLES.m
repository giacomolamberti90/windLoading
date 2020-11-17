function tile = load_charLES(angle, label, path)

% taps list - reference
tapsList   = load('/home/giacomol/Desktop/Research/windLoading/windTunnel/matlabCode/ArupGUI3.0/Orifizi-v2.mat');
taps_A_ref = cell2mat(tapsList.nomi{1}');
taps_B_ref = cell2mat(tapsList.nomi{2}');

% remove corrupted data
[~, ~, index] = setxor({'B1315'}, taps_B_ref);
taps_B_ref = taps_B_ref(index,:);

%% load data
pbin = fullfile(path, 'test.pbin'); % xyz data 
fh = fopen(pbin, 'r');
endian = 'l';

magic_number = fread(fh, 1, 'int64', endian); % magic number
if magic_number ~= 1235813 
    endian = 'b';
end

fread(fh, 1, 'int64', endian); % skip version
np = fread(fh, 1,'int64', endian); % number of points
assert(2 == fread(fh, 1, 'int64', endian)); % 0:no data, 1:delta, 2:index
buf = fread(fh, 3*np, 'double', endian); % always double for consistency w/ other pbins
ind = fread(fh, np, 'int64', endian); % global index (should equal ordering in ascii file)
fclose(fh);

xyz = zeros(np, 3);
for i = 1:np
    for j = 1:3
        xyz(1 + ind(i), j) = buf(3*(i-1) + j);
    end
end

theta = -angle * pi/180;
for i = 1:np
    for j = 1:3
        roxyz(i,1) = (xyz(i,1) - 0.5) * cos(theta) + xyz(i,3) * sin(theta) + 0.5;  % 45
        roxyz(i,3) = (xyz(i,1) - 0.5) * sin(theta) - xyz(i,3) * cos(theta);        % 45
        roxyz(i,2) = xyz(i,2);
    end
end

t = linspace(0, 100, 200000);
pvalue = zeros(200000, np);
idx = 0;

for irk = 1:200000
    
    tag = num2str(irk, '%08d');
    pcd = strcat(path, 'test.', tag, '.pcd'); % scalar data    

    fh = fopen(pcd, 'r');
    
    if fh ~= -1
        idx = idx + 1;
        
        endian = 'l';
        magic_number = fread(fh, 1, 'int64', endian); % magic number

        if magic_number ~= 1235813 
            endian = 'b';
        end

        fread(fh,1,'int64',endian); % skip version
        assert(np == fread(fh,1,'int64',endian)); % number of points
        nv = fread(fh,1,'int64',endian); % number of variables
        prec = fread(fh,1,'int64',endian); % 0 float/ 1 double

        if (prec == 1)
            buf = fread(fh,nv*np,'double',endian); % note that each scalar gets a row
        else
            buf = fread(fh,nv*np,'float',endian);
        end

        scd = zeros(nv,np);
        for j = 1:nv
            for i = 1:np
                scd(j,1+ind(i)) = buf(i+(j-1)*np);
            end
        end
        fclose(fh);
        pvalue(idx, :) = scd(nv,:)/36.7;
    end
end

%% define data structure
tile = {};

% tileA
if label == 'A'
    tile.CI          = 'no';
    tile.taps        = taps_A_ref(1:end-1,:);
    tile.coords      = roxyz(:, 1:2);
    tile.time        = t;
    tile.timeHistory = pvalue;
    tile.mean        = mean(pvalue)';
    tile.std         = std(pvalue)';
    tile.U           = 7.74;
    tile.q           = 36.7;
elseif  label == 'B'
    tile.CI          = 'no';
    tile.taps        = taps_B_ref;
    tile.coords      = roxyz(:, 1:2);
    tile.time        = t;
    tile.timeHistory = pvalue;
    tile.mean        = mean(pvalue)';
    tile.std         = std(pvalue)';
    tile.U           = 7.74;
    tile.q           = 36.7;
end
