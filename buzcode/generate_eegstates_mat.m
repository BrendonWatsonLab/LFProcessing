function generate_eegstates_mat(basePath, channels)
% generate_eegstates_mat - Generate baseName.eegstates.mat file containing spectrogram
%
%  USAGE
%    generate_eegstates_mat(basePath, channels)
%
%    basePath  - path to directory containing baseName.xml and baseName.eeg/.lfp
%    channels  - vector of channel numbers (0-based indexing) to include in spectrogram.
%
%  NOTES:
%    - Requires LoadXml.m and AR/spectrogram helpers from original code.
%    - Assumes .eeg or .lfp file has same base name as .xml
%    - Adjust spectrogram parameters (WinLength, NW, etc.) as needed.
%
%  OUTPUT:
%    Creates baseName.eegstates.mat with a StateInfo struct:
%       StateInfo.Chs       = chosen channels
%       StateInfo.nCh       = number of channels chosen
%       StateInfo.motion    = [] or motion signal if available
%       StateInfo.eegFS     = sampling frequency (Hz)
%       StateInfo.fspec     = cell array of spectrogram data structures
%
%  The spectrogram data structure fspec{i} includes fields:
%       fspec{i}.spec - power spectra
%       fspec{i}.fo   - frequency axis
%       fspec{i}.to   - time axis
%
%  After running, you can later load the .eegstates.mat in StateEditor.


%% Identify base name and files
[basePath, baseName, ~] = fileparts(basePath);
if isempty(basePath), basePath = pwd; end
xmlFile = fullfile(basePath, [baseName '.xml']);
eegFile = fullfile(basePath, [baseName '.eeg']);
if ~exist(eegFile,'file')
    eegFile = fullfile(basePath, [baseName '.lfp']);
    if ~exist(eegFile,'file')
        error('No .eeg or .lfp file found.');
    end
end

%% Load metadata from XML
if ~exist(xmlFile,'file')
    error('No corresponding .xml file found.');
end
Par = LoadXml(xmlFile);
nChannels = Par.nChannels;
eegFS = Par.rates.wideband; % or Par.SampleRate if defined

% Validate chosen channels
if any(channels<0) || any(channels>=nChannels)
    error('Chosen channels out of range. Must be between 0 and %d.', nChannels-1);
end

%% Load EEG data from chosen channels
disp('Loading EEG/LFP data...');
fid = fopen(eegFile,'r');
if fid == -1
    error('Unable to open eeg/lfp file.');
end
fseek(fid,0,'eof');
fSize = ftell(fid);
fseek(fid,0,'bof');
samples = fSize/(nChannels*2); % int16=2 bytes
samples = floor(samples);
dataMatrix = fread(fid,[nChannels samples],'int16');
fclose(fid);
dataMatrix = double(dataMatrix(channels+1,:)); % +1 since channels are 0-based
% Now dataMatrix is [length(channels) x samples]

%% Whiten and compute spectrogram
disp('Whitening and computing spectrogram...');
% Whiten each channel
weeg = cell(1,length(channels));
for i = 1:length(channels)
    weeg{i} = WhitenSignalIn(dataMatrix(i,:)',eegFS*2000,1); % adjust params as needed
end

% Convert cell of whitened channels into a matrix
wmat = cell2mat(weeg'); % [length(channels) x samples]
wmat = wmat'; % make it samples x channels for mtchglongIn if needed

% Choose spectrogram params
WinLength = 3072; % example
NW = 3; 
nFFT = WinLength;
nOverlap = 0;   % no overlap for simplicity, adjust as needed
FreqRange = [0 200]; % freq range

% Compute spectrogram for each channel
fspec = cell(1,length(channels));
for iCh = 1:length(channels)
    [spec,fo,to] = mtchglongIn(weeg{iCh}, nFFT, eegFS, WinLength, nOverlap, NW, [], [], FreqRange);
    fspec{iCh}.spec = single(spec);
    fspec{iCh}.fo   = fo;
    fspec{iCh}.to   = to;
    fspec{iCh}.info.Ch = channels(iCh);
    fspec{iCh}.info.FileInfo.name = eegFile;
end

%% If you have motion data, load/compute it here
% For now, we skip motion or set to empty
motion = [];

%% Construct StateInfo struct
disp('Constructing StateInfo and saving...');
StateInfo = struct();
StateInfo.nCh = length(channels);
StateInfo.Chs = channels;
StateInfo.eegFS = eegFS;
StateInfo.motion = motion;
StateInfo.fspec = fspec;

% Optionally, save rawEeg in StateInfo for portability
% This will make the file large. If you want:
% StateInfo.rawEeg = dataMatrix(channels+1,:); % original data if desired

save(fullfile(basePath, [baseName '.eegstates.mat']), 'StateInfo');
disp('Done. baseName.eegstates.mat created.');
end
