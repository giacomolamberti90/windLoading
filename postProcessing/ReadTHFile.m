function varargout = ReadTHFile(varargin)

%ReadTHFile   Reads a time-history data file
%
%  For VELOCITY DATA time-history files:
%    [U, V, W, Ps, Settings] = ReadTHFile(fileName, numSamplesToProcess, offset);
%
%    Returns U, V, W, Pstatic and settings
%    Also returns reference pressure data if sampled
%
%
%  For ANALOGUE DATA time-history files:
%    [Data, Settings] = ReadTHFile(fileName, numSamplesToProcess, offset);
%
%    Returns data (one column per channel) and settings
%
%
%  Note: File name is optional - enter '' to get a file open dialog box
%	      Number of samples to process is optional - enter '0' to read the entire file
%        Offset is optional - enter number of samples to offset by
%
%See also VELPITCHYAW, TRANSFORMAXES and ALIGNWITHFLOW


tempFileName = '';
numSamplesToProcess = 0;
offset = 0;

% Get input arguments
if (nargin >= 1)
   tempFileName = varargin{1};
end
if (nargin >= 2)
   numSamplesToProcess = varargin{2}; 
end
if (nargin >= 3)
   offset = varargin{3};
end

% Show file open dialog if required
if (length(tempFileName) == 0)
   [filename,filepath]=uigetfile('*.th?;*.th;*.bp?','Select time-history file to open...');
   if (filename == 0)
     fprintf('No file selected!\n')
     return
   end;
   if (filename == 0)
     fprintf('No path selected!\n')
     return
   end;
   
   file=filepath;
   file((length(filepath)+1):(length(filepath)+length(filename))) = filename;
else
   file = tempFileName;
end;

% Open time-history file
fid = fopen(file, 'r');

if fid >= 3
   %	Read header
   fileFormat = fread(fid,1,'int32');
%   fprintf('File format : %i\n', fileFormat)
   
   if (fileFormat == 1) | (fileFormat == 6401)
      probeID = fread(fid,1,'int32');
      deviceName = sprintf('Probe %0.3i', deviceID);
      dataType = 100;
      dateString = '';
      timeString = '';
      numSamples = fread(fid,1,'int32');
      blockSize = fread(fid,1,'int32');
      dataRate = fread(fid,1,'float64');
      hasPref = fread(fid,1,'uint8');
   
      if (numSamples == 0)
         numSamples = NumSamplesPerChannel(fid,numChannels);
      end
   
      fprintf('Probe ID : %0.3i\n', probeID)
      fprintf('Num samples : %i\n', numSamples)
      fprintf('Data rate : %0.1f Hz\n', dataRate)
      if (hasPref)
         fprintf('Has Pref : yes\n')
      else
         fprintf('Has Pref : no\n')
      end;
   elseif fileFormat == 2
      probeID = fread(fid,1,'int32');
      deviceName = sprintf('Probe %0.3i', deviceID);
      dataType = 100
      dateString = '';
      timeString = '';
      numSamples = fread(fid,1,'int32');
      blockSize = fread(fid,1,'int32');
      dataRate = fread(fid,1,'float64');
      Pbaro = fread(fid,1,'float64');
      meanTemperature = fread(fid,1,'float64');
      hasPref = fread(fid,1,'uint8');
      
      if (numSamples == 0)
         numSamples = NumSamplesPerChannel(fid,numChannels);
      end
   
      fprintf('Probe ID : %0.3i\n', probeID)
      fprintf('Num samples : %i\n', numSamples)
      fprintf('Data rate : %0.1f Hz\n', dataRate)
      fprintf('Barometric pressure : %0.1f Pa\n', Pbaro)
      fprintf('Mean temperature : %0.1f °C\n', meanTemperature)
      if (hasPref)
         fprintf('Has ref pressure : yes\n')
      else
         fprintf('Has ref pressure : no\n')
      end;
   elseif fileFormat == 3
      deviceType = fread(fid,1,'int32');
      deviceID = fread(fid,1,'int32');
      dataType = fread(fid,1,'int32');
      deviceName = sprintf('%s %0.3i', DeviceTypeName(deviceType), deviceID);
      dateString = ReadDate(fid);
      timeString = ReadTime(fid);
      numSamples = fread(fid,1,'int32');
      blockSize = fread(fid,1,'int32');
      dataRate = fread(fid,1,'float64');
      Pbaro = fread(fid,1,'float64');
      meanTemperature = fread(fid,1,'float64');
      hasPref = fread(fid,1,'uint8');
      
      if (numSamples == 0)
         numSamples = NumSamplesPerChannel(fid,numChannels);
      end
   
      fprintf('Device name         : %s\n', deviceName)
      fprintf('Num samples         : %i\n', numSamples)
      fprintf('Data rate           : %0.1f Hz\n', dataRate)
      fprintf('Barometric pressure : %0.1f Pa\n', Pbaro)
      fprintf('Mean temperature    : %0.1f °C\n', meanTemperature)
      if (hasPref)
         fprintf('Has ref pressure    : yes\n')
      else
         fprintf('Has ref pressure    : no\n')
      end;
   elseif fileFormat == 101
      dataType = fread(fid,1,'int32');
      dataOrder = fread(fid,1,'int32');
      numChannels = fread(fid,1,'int32');
      numSamples = fread(fid,1,'int32');
      blockSize = fread(fid,1,'int32');
      dataRate = fread(fid,1,'float64');
      dateString = ReadDate(fid);
      timeString = ReadTime(fid);
      for n = 1:numChannels; chanLabels(n) = {sprintf('Ch %d',n)}; end
      for n = 1:numChannels; chanUnitAbrs(n) = {''}; end
      
      if (numSamples == 0)
         numSamples = NumSamplesPerChannel(fid,numChannels);
      end
   
      fprintf('Num channels : %i\n', numChannels)
      fprintf('Num samples  : %i\n', numSamples)
      fprintf('Data rate    : %0.1f Hz\n', dataRate)
   elseif fileFormat == 102
      headerSize = fread(fid,1,'int32');
      dataType = fread(fid,1,'int32');
      dataOrder = fread(fid,1,'int32');
      numChannels = fread(fid,1,'int32');
      numSamples = fread(fid,1,'int32');
      blockSize = fread(fid,1,'int32');
      dataRate = fread(fid,1,'float64');
      dateString = ReadDate(fid);
      timeString = ReadTime(fid);
      for n = 1:numChannels
         labelLength = fread(fid,1,'int32');
         chanLabels(n) = { char( fread(fid,labelLength,'uchar')' ) };
      end
      for n = 1:numChannels
         unitAbrLength = fread(fid,1,'int32');
         chanUnitAbrs(n) = { char( fread(fid,unitAbrLength,'uchar')' ) };
      end
      
      if (numSamples == 0)
         numSamples = NumSamplesPerChannel(fid,numChannels);
      end
   
      fprintf('Num channels : %i\n', numChannels)
      fprintf('Num samples  : %i\n', numSamples)
      fprintf('Data rate    : %0.1f Hz\n', dataRate)
   elseif fileFormat == 201
      deviceType = fread(fid,1,'int32');
      deviceID = fread(fid,1,'int32');
      deviceName = sprintf('%s %0.3i', DeviceTypeShortName(deviceType), deviceID);
      dataType = fread(fid,1,'int32');
      dataOrder = fread(fid,1,'int32');
      numChannels = fread(fid,1,'int32');
      numSamples = fread(fid,1,'int32');
      blockSize = fread(fid,1,'int32');
      dataRate = fread(fid,1,'float64');
      dateString = ReadDate(fid);
      timeString = ReadTime(fid);
      for n = 1:numChannels; chanLabels(n) = {sprintf('Ch %d',n)}; end
      for n = 1:numChannels; chanUnitAbrs(n) = {''}; end
      
      if (numSamples == 0)
         numSamples = NumSamplesPerChannel(fid,numChannels);
      end
   
      fprintf('Device name  : %s\n', deviceName)
      fprintf('Num channels : %i\n', numChannels)
      fprintf('Num samples  : %i\n', numSamples)
      fprintf('Data rate    : %0.1f Hz\n', dataRate)
   elseif fileFormat == 202
      headerSize = fread(fid,1,'int32');
      deviceNameLength = fread(fid,1,'int32');
      deviceName = char( fread(fid,deviceNameLength,'uchar')' );
      deviceType = fread(fid,1,'int32');
      deviceID = fread(fid,1,'int32');
      dataType = fread(fid,1,'int32');
      dataOrder = fread(fid,1,'int32');
      numChannels = fread(fid,1,'int32');
      numSamples = fread(fid,1,'int32');
      blockSize = fread(fid,1,'int32');
      dataRate = fread(fid,1,'float64');
      dateString = ReadDate(fid);
      timeString = ReadTime(fid);
      for n = 1:numChannels
         labelLength = fread(fid,1,'int32');
         chanLabels(n) = { char( fread(fid,labelLength,'uchar')' ) };
      end
      for n = 1:numChannels
         unitAbrLength = fread(fid,1,'int32');
         chanUnitAbrs(n) = { char( fread(fid,unitAbrLength,'uchar')' ) };
      end
      
      if (numSamples == 0)
         numSamples = NumSamplesPerChannel(fid,numChannels);
      end
   
      fprintf('Device name  : %s\n', deviceName)
      fprintf('Num channels : %i\n', numChannels)
      fprintf('Num samples  : %i\n', numSamples)
      fprintf('Data rate    : %0.1f Hz\n', dataRate)
   else
      disp('File format not supported!	')
      return      
   end
   
   % Configure settings output
   if (fileFormat == 1) | (fileFormat == 6401) | (fileFormat == 2) | (fileFormat == 3)
      varargout{5} = struct('DeviceName',deviceName, 'DeviceType',DeviceTypeName(deviceType), 'DeviceID',deviceID, ...
                            'DataType',DataTypeName(dataType), ...
                            'BlockSize',blockSize, 'DataRate',dataRate, ...
                            'Pbaro',Pbaro, 'Tmean',meanTemperature, ...
                            'Pref',[], 'FirstSampleDate',dateString, 'FirstSampleTime',timeString);
   elseif (fileFormat == 101) | (fileFormat == 102)
      varargout{2} = struct('DataType',DataTypeName(dataType), 'DataOrder',DataOrderName(dataOrder), ...
                            'NumChannels',numChannels, 'NumSamples',numSamples, 'BlockSize',blockSize, 'DataRate',dataRate, ...
                            'FirstSampleDate',dateString, 'FirstSampleTime',timeString, ...
                            'ChanLabels',{chanLabels}', 'ChanUnitAbrs',{chanUnitAbrs});
   elseif (fileFormat == 201) | (fileFormat == 202)
      varargout{2} = struct('DeviceName',deviceName, 'DeviceType',DeviceTypeName(deviceType), 'DeviceID',deviceID, ...
                            'DataType',DataTypeName(dataType), 'DataOrder',DataOrderName(dataOrder), ...
                            'NumChannels',numChannels, 'NumSamples',numSamples, 'BlockSize',blockSize, 'DataRate',dataRate, ...
                            'FirstSampleDate',dateString, 'FirstSampleTime',timeString, ...
                            'ChanLabels',{chanLabels}', 'ChanUnitAbrs',{chanUnitAbrs});
   end
   
   fprintf('\n');
      
   
   % Read data
   if (fileFormat == 1) | (fileFormat == 6401) | (fileFormat == 2) | (fileFormat == 3)
      % Seek to required start position
      if (offset > numSamples)
         disp('Offset larger than number of samples in file!\n')
         return
      end
      blockStartSample = floor(offset/blockSize)*blockSize;
      if (hasPref)
         fseek(fid, blockStartSample*4*5, 'cof');
      else
         fseek(fid, blockStartSample*4*4, 'cof');
      end
      
      if (numSamplesToProcess == 0)
         numBlocks = numSamples / blockSize;
      else
         numBlocks = floor(numSamplesToProcess/blockSize);
      end;

      varargout{1} = zeros(numBlocks*blockSize,1);
      varargout{2} = zeros(numBlocks*blockSize,1);
      varargout{3} = zeros(numBlocks*blockSize,1);
      varargout{4} = zeros(numBlocks*blockSize,1);
      if (hasPref)
         varargout{5}.Pref = zeros(numBlocks*blockSize,1);
      end;
   
      for blockNum = 1:numBlocks
         index1 = (blockNum-1) * blockSize + 1;
         index2 = blockNum * blockSize;
      
         varargout{1}(index1:index2) = fread(fid,blockSize,'float32');
         varargout{2}(index1:index2) = fread(fid,blockSize,'float32');
         varargout{3}(index1:index2) = fread(fid,blockSize,'float32');
         varargout{4}(index1:index2) = fread(fid,blockSize,'float32');
         if (hasPref)
            varargout{5}.Pref(index1:index2) = fread(fid,blockSize,'float32');
         end;
      end;
   
   elseif (fileFormat == 101) | (fileFormat == 102) | (fileFormat == 201) | (fileFormat == 202)
      % Seek to required start position
      if (offset > numSamples)
         disp('Offset larger than number of samples in file!\n')
         return
      end
      blockStartSample = floor(offset/blockSize)*blockSize;
      fseek(fid, blockStartSample*numChannels*4, 'cof');
      
      if (numSamplesToProcess == 0)
         numBlocks = numSamples / blockSize;
      else
         numBlocks = floor(numSamplesToProcess/blockSize);
      end;

      if dataOrder == 0			% doNonInterleaved
         varargout{1} = zeros(numBlocks*blockSize, numChannels);
         for blockNum = 1:numBlocks
            index1 = (blockNum-1) * blockSize + 1;
            index2 = blockNum * blockSize;
      
            varargout{1}(index1:index2,:) = fread(fid,[blockSize,numChannels],'float32');
         end;
      elseif dataOrder == 1	% doInterleaved
         varargout{1} = fread(fid,[numChannels,blockSize],'float32')';
      end;
   
   end;   
   

   % Close file
   fclose(fid);
else
   fprintf('Failed to open file!\n\n');
end;

return



% Convert device type to device type name
function name = DeviceTypeName(DeviceType)

switch DeviceType
   case 0,
      name = 'Cobra';
   case 4,
      name = 'Four-hole Cobra';
   case 5,
      name = 'Five-hole Cobra';
   case 13,
      name = 'Thirteen-hole ECA';
   case 1000,
      name = 'DP Module';
   case 1015,
      name = '15-channel DP Module';
   case 1016,
      name = '16-channel DP Module';
   case 1032,
      name = '32-channel DP Module';
   case 1064,
      name = '64-channel DP Module';
   case 1128,
      name = '128-channel DP Module';
   case 1256,
      name = '256-channel DP Module';
   case 8000,
      name = 'Force Balance';
   case 8010,
      name = 'JR3 Balance';
   case 8011.
      name = 'JR3 Balance (A/D)';
   case 8012,
      name = 'JR3 Balance (DSP)';
   case 8020,
      name = 'Aeroelastic Base Balance';
   case 10000,
      name = 'Analog Data Device';
   otherwise
      name = 'Unknown';
end
   
return



% Convert device type to device type short name
function name = DeviceTypeShortName(DeviceType)

switch DeviceType
   case {0, 4, 5},
      name = 'Cobra';
   case 13,
      name = 'ECA';
   case {1000, 1015, 1016, 1032, 1064, 1128, 1256},
      name = 'DPM';
   case 8000,
      name = 'Force Balance';
   case {8010, 8011, 8012},
      name = 'JR3';
   case 8020,
      name = 'ABB';
   case 10000,
      name = 'ADD';
   otherwise
      name = 'Unknown';
end
   
return



% Convert data type to data type name
function name = DataTypeName(DataType)

switch DataType
   case 0,
      name = 'Analogue';
   case 10,
      name = 'Voltage';
   case 20,
      name = 'Pressure';
   case 30,
      name = 'Temperature';
   case 100,
      name = 'Velocity';
   case 200,
      name = 'Force-Moment';
   otherwise
      name = 'Unknown';
end
   
return



% Convert data order to data order name
function name = DataOrderName(DataOrder)

switch DataOrder
   case 0,
      name = 'Non-interleaved';
   case 1,
      name = 'Interleaved';
   otherwise
      name = 'Unknown';
end
   
return



% Read a date and time
function dateTimeString = ReadDateTime(fid)

dateString = ReadDate(fid);
timeString = ReadTime(fid);

dateTimeString = sprintf('%s  %s', dateString, timeString);

return
   
   
   
% Read a date
function dateString = ReadDate(fid)

year = fread(fid, 1, 'int16');
month = fread(fid, 1, 'int16');
dayOfWeek = fread(fid, 1, 'int16');
day = fread(fid, 1, 'int16');

dateNumber = datenum(year, month, day);
dateString = datestr(dateNumber, 'dd-mmm-yyyy');

return
   
   
   
% Read a time
function timeString = ReadTime(fid)

hour = fread(fid, 1, 'int16');
minute = fread(fid, 1, 'int16');
second = fread(fid, 1, 'int16');
milliSecond = fread(fid, 1, 'int16');

timeNumber = datenum(0, 0, 0, hour, minute, second);
timeString = datestr(timeNumber, 'HH:MM:SS');
msString = sprintf('%0.3f', milliSecond/1000);
msString = msString(2:end);
timeString = strcat(timeString, msString);

return



% Determine number of samples per channel
function numSamples = NumSamplesPerChannel(fid, numChannels)

   filePos = ftell(fid);
   fseek(fid,0,'eof');
   fileLength = ftell(fid);
   fseek(fid,filePos,'bof');
   numSamples = (fileLength - filePos) / (numChannels*4);

return