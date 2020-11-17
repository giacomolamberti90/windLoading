function [RawData]=importScanivalveDataFile(fName)    
% Reads binary data files from pressure scanning system
% fName: Text string containing data file name and path 
% RawData: Matrix of scanned differential pressures in psi. Each column 
% represents a pressure tap. Each row represents a scanned frame. 
% Note: Saturated transducers output 9999 psi. This code converts them to
% NaN.
FID = fopen(fName);
raw = fread(FID,'single');
fclose(FID);
raw = single(raw);
% read the number of channels and determine the number of samples
nChannels = (typecast(raw(1),'uint16'));
nChannels = double(nChannels(2));
nSamples = double(length(raw)/(1+(nChannels+2)));
% parse the data, so that each column belongs to a channel
ss = 1:nSamples;
RawData=zeros(nSamples,nChannels);
for cc = 1:nChannels
    RawData(ss,cc) = double(raw(3+(ss-1)*(nChannels+3)+cc));
end
% remove over voltages
[a,b] = find(RawData == 9999);
for jj = 1:length(a)
    RawData(a(jj),b(jj)) = nan;
end
    
        
        