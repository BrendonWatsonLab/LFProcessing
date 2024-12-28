function splitDatFile(datFilePath, outputDir)
% splitDatFile - Split a large .dat file into separate files, one per channel.
%
%   This function reads a .dat file and its associated .xml metadata file.
%   It then splits the data into individual channels and saves each channel's
%   data as a separate binary (.bin) or MATLAB (.mat) file.
%
%   Inputs:
%       datFilePath - Path to the .dat file.
%       outputDir   - Directory where individual channel files will be saved.
%
%   Example:
%       splitDatFile('data.dat', 'outputChannels');
%
%   Author: Your Name
%   Date: YYYY-MM-DD

    %% Parse File and Metadata
    [dataDir, baseName, ext] = fileparts(datFilePath);
    if ~strcmp(ext, '.dat')
        error('Provided file must be a .dat file.');
    end
    xmlFilePath = fullfile(dataDir, [baseName '.xml']);
    if ~exist(xmlFilePath, 'file')
        error('The corresponding .xml file does not exist: %s', xmlFilePath);
    end

    % Extract metadata from .xml
    doc = xmlread(xmlFilePath);
    nChannels = str2double(char(doc.getElementsByTagName('nChannels').item(0).getFirstChild.getData));
    samplingRate = str2double(char(doc.getElementsByTagName('samplingRate').item(0).getFirstChild.getData));

    % Validate output directory
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    %% Memory Map the File
    fileInfo = dir(datFilePath);
    totalSamples = fileInfo.bytes / (2 * nChannels); % int16 data, 2 bytes per sample

    m = memmapfile(datFilePath, 'Format', 'int16');

    %% Split Data into Individual Channel Files
    fprintf('Splitting %d channels...\n', nChannels);
    for ch = 1:nChannels
        % Extract data for the current channel
        channelData = double(m.Data(ch:nChannels:end));

        % Option 1: Save as binary file
        binFileName = fullfile(outputDir, sprintf('%s_channel%d.bin', baseName, ch));
        fileID = fopen(binFileName, 'w');
        fwrite(fileID, channelData, 'double');
        fclose(fileID);

        % Option 2: Save as MATLAB .mat file (optional, larger files)
        % matFileName = fullfile(outputDir, sprintf('%s_channel%d.mat', baseName, ch));
        % save(matFileName, 'channelData', '-v7.3'); % Use -v7.3 for large arrays

        fprintf('Channel %d saved to %s\n', ch, binFileName);
    end

    fprintf('Splitting complete. Data saved in %s.\n', outputDir);
end
