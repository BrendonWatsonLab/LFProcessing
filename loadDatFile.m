function [dataMatrix, timeVector, metaData] = loadDatFile(datFilePath)
% loadDatFile - Load raw LFP data from a .dat file recorded with Neuroscope.
%
%   This function reads the specified .dat file and its corresponding .xml
%   metadata file (expected to have the same base name and reside in the
%   same directory). It then interprets the raw, interleaved data according
%   to the parameters found in the .xml file, separating it into individual
%   channels. Unlike previous versions, if the total number of data points 
%   is not divisible by the number of channels specified in the .xml file, 
%   it trims the excess data to match the channel count.
%
%   Outputs:
%       dataMatrix: A [nChannels x nSamples] double array containing raw data.
%                   Each row corresponds to one channel.
%       timeVector: A [1 x nSamples] double vector of timestamps (seconds).
%       metaData:   A struct with metadata fields:
%                     - filePath
%                     - xmlPath
%                     - numChannels
%                     - samplingRate
%                     - numSamples
%                     - recordingTimeSec
%
%   Requirements:
%       * The .xml file must contain <nChannels> and <samplingRate> tags.
%       * If total samples do not perfectly divide into channels, the extra 
%         samples at the end will be discarded.
%
%   Example:
%       [dataMatrix, timeVector, metaData] = loadDatFile('path/to/yourData.dat');
%
%   Author: Your Name
%   Date: YYYY-MM-DD
%   -------------------------------------------------------------

    %% -------------------- Parameter & File Handling --------------------
    [dataDir, baseName, ext] = fileparts(datFilePath);
    if ~strcmp(ext, '.dat')
        error('Provided file must be a .dat file.');
    end
    xmlFilePath = fullfile(dataDir, [baseName '.xml']);
    
    if ~exist(datFilePath, 'file')
        error('The specified .dat file does not exist: %s', datFilePath);
    end
    if ~exist(xmlFilePath, 'file')
        error('The corresponding .xml file does not exist: %s', xmlFilePath);
    end

    %% -------------------- Parse the XML File --------------------
    doc = xmlread(xmlFilePath);

    % Extract nChannels
    nChannelsNode = doc.getElementsByTagName('nChannels');
    if nChannelsNode.getLength == 0
        error('No nChannels tag found in XML.');
    end
    nChannels = str2double(char(nChannelsNode.item(0).getFirstChild.getData));

    % Extract samplingRate
    sampleRateNode = doc.getElementsByTagName('samplingRate');
    if sampleRateNode.getLength == 0
        error('No samplingRate tag found in XML.');
    end
    samplingRate = str2double(char(sampleRateNode.item(0).getFirstChild.getData));

    %% -------------------- Load Binary Data from the .dat File --------------------
    fileID = fopen(datFilePath, 'r');
    if fileID == -1
        error('Unable to open .dat file: %s', datFilePath);
    end

    rawData = fread(fileID, 'int16');
    fclose(fileID);

    totalSamples = length(rawData);
    remainder = mod(totalSamples, nChannels);
    if remainder ~= 0
        % Trim the excess data
        truncatedSamples = totalSamples - remainder;
        warning(['Total data points (', num2str(totalSamples), ') not divisible by channels (', ...
                 num2str(nChannels), '). Trimming the last ', num2str(remainder), ' samples.']);
        rawData = rawData(1:truncatedSamples);
        totalSamples = truncatedSamples;
    end

    nSamples = totalSamples / nChannels;

    dataMatrix = reshape(rawData, [nChannels, nSamples]);
    dataMatrix = double(dataMatrix);  % Convert to double precision

    %% -------------------- Create Time Vector --------------------
    timeVector = (0 : nSamples-1) / samplingRate;

    %% -------------------- Construct MetaData Struct --------------------
    metaData = struct();
    metaData.filePath         = datFilePath;
    metaData.xmlPath          = xmlFilePath;
    metaData.numChannels      = nChannels;
    metaData.samplingRate     = samplingRate;
    metaData.numSamples       = nSamples;
    metaData.recordingTimeSec = nSamples / samplingRate;

end


