function [dataMatrix, timeVector, metaData] = loadFilteredLFP(datFilePath, outFs, lopass)
% loadFilteredLFP - Load a .dat file and produce a filtered, downsampled LFP-like dataset.
%
%   This function is inspired by the lab code (bz_LFPfromDat) and 
%   the original loadDatFile function. It:
%       - Reads the raw .dat file in chunks.
%       - Applies a low-pass filter (using a sinc filter) to extract LFP range.
%       - Downsamples the data from the original sampling rate (inFs) to a lower rate (outFs).
%       - Returns the resulting data in memory as dataMatrix, along with timeVector and metaData.
%
%   Inputs:
%       datFilePath - Path to the .dat file
%       outFs       - (optional) Desired output sampling rate of the LFP. Default: 1250 Hz
%       lopass      - (optional) Low-pass cutoff frequency. Default: 450 Hz
%
%   Outputs:
%       dataMatrix: [nChannels x nSamples_downsampled] double matrix of LFP data
%       timeVector: [1 x nSamples_downsampled] time vector in seconds at outFs
%       metaData:   struct with fields:
%                      * filePath
%                      * xmlPath
%                      * numChannels
%                      * originalFs   (inFs)
%                      * samplingRate (outFs)
%                      * numSamples   (after downsampling)
%                      * recordingTimeSec (after downsampling)
%                      * lowPassFreq  (lopass)
%
%   Requirements:
%       - The .xml file (same basename) to extract nChannels and samplingRate.
%       - iosr toolbox for iosr.dsp.sincFilter (https://github.com/IoSR-Surrey/MatlabToolbox)
%
%   Example:
%       [dataMatrix, timeVector, metaData] = loadFilteredLFP('mydata.dat', 1250, 450);
%
%   Note:
%       This code processes the entire recording into memory. For extremely large files, consider 
%       modifying to process segments or saving intermediate results. Also, ensure you have enough RAM.
%
%   Author: Your Name
%   Date: YYYY-MM-DD
%   -------------------------------------------------------------

    if ~exist('outFs','var') || isempty(outFs)
        outFs = 1250;  % default LFP sampling rate
    end
    if ~exist('lopass','var') || isempty(lopass)
        lopass = 450;  % default low pass cutoff
    end

    % Check file
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

    %% Parse XML for metadata
    doc = xmlread(xmlFilePath);
    nChannelsNode = doc.getElementsByTagName('nChannels');
    if nChannelsNode.getLength == 0
        error('No nChannels tag found in XML.');
    end
    nChannels = str2double(char(nChannelsNode.item(0).getFirstChild.getData));

    sampleRateNode = doc.getElementsByTagName('samplingRate');
    if sampleRateNode.getLength == 0
        error('No samplingRate tag found in XML.');
    end
    inFs = str2double(char(sampleRateNode.item(0).getFirstChild.getData));

    % Check Nyquist
    if lopass > inFs/2
        warning('Low pass cutoff beyond Nyquist of original data. Adjusting lopass.');
        lopass = inFs/2;
    end

    % Prepare to read .dat file
    fileInfo = dir(datFilePath);
    nBytes = fileInfo.bytes;
    sizeInBytes = 2; % int16 = 2 bytes
    totalSamples = nBytes / (sizeInBytes * nChannels);
    totalTime = totalSamples / inFs;

    % Parameters for filtering and downsampling
    % ratio = lopass/(inFs/2)
    ratio = lopass / (inFs/2);
    sampleRatio = inFs / outFs;

    if mod(sampleRatio,1)~=0
        warning('inFs/outFs is not an integer. Consider choosing outFs that divides inFs.');
    end

    % Decide on chunk sizes
    % Use a chunk size that's divisible by sampleRatio for convenience.
    chunksize = 1e5; 
    if mod(chunksize, sampleRatio) ~= 0
        chunksize = chunksize + sampleRatio - mod(chunksize, sampleRatio);
    end

    % ntbuff: buffer for filter
    ntbuff = 525; 
    if mod(ntbuff, sampleRatio)~=0
        ntbuff = ntbuff + sampleRatio - mod(ntbuff, sampleRatio);
    end

    % Number of full chunks
    nbChunks = floor(totalSamples/chunksize);

    % Open the file
    fidI = fopen(datFilePath, 'r');
    if fidI == -1
        error('Unable to open .dat file.');
    end

    % We will store all processed data in memory
    % Estimate final number of samples after downsampling:
    % Each chunk of size 'chunksize' will produce about (chunksize/sampleRatio) samples after downsampling
    % There's a buffer effect at the start, but approximately:
    estSamplesPerChunk = chunksize/sampleRatio;
    estTotalSamplesDown = nbChunks*estSamplesPerChunk; % plus remainder
    % We'll store data in a growing array (or preallocate if memory allows)
    % To be safe, let's first store in a cell and concatenate later if huge.

    % Preallocate large arrays might be huge. We'll just store in a cell and cat at end.
    dataSegments = cell(nbChunks+1,1);

    import iosr.dsp.* % Ensure iosr toolbox is on path

    % Process each chunk
    for ibatch = 1:nbChunks
        if ibatch > 1
            fseek(fidI, ((ibatch-1)*nChannels*sizeInBytes*chunksize)-(nChannels*sizeInBytes*ntbuff), 'bof');
            dat = fread(fidI, nChannels*(chunksize+2*ntbuff), 'int16');
            dat = reshape(dat, [nChannels (chunksize+2*ntbuff)]);
        else
            % First batch
            dat = fread(fidI, nChannels*(chunksize+ntbuff), 'int16');
            dat = reshape(dat, [nChannels (chunksize+ntbuff)]);
        end

        downData = nan(nChannels, chunksize/sampleRatio);
        for ch = 1:nChannels
            d = double(dat(ch,:));
            tmp = sincFilter(d, ratio);
            if ibatch == 1
                % first batch: use samples from sampleRatio:sampleRatio:end-ntbuff
                validSamples = tmp(sampleRatio:sampleRatio:end-ntbuff);
            else
                % subsequent batches: skip initial ntbuff samples
                validSamples = tmp(ntbuff+sampleRatio:sampleRatio:end-ntbuff);
            end
            downData(ch,:) = validSamples;
        end

        dataSegments{ibatch} = downData;
    end

    % Remainder
    remainder = totalSamples - nbChunks*chunksize;
    if remainder > 0
        fseek(fidI, (nbChunks*chunksize*nChannels*sizeInBytes)-(nChannels*sizeInBytes*ntbuff), 'bof');
        dat = fread(fidI, nChannels*(remainder+ntbuff), 'int16');
        dat = reshape(dat,[nChannels (remainder+ntbuff)]);

        remainderSamples = floor(remainder/sampleRatio);
        if remainderSamples > 0
            downData = nan(nChannels, remainderSamples);
            for ch = 1:nChannels
                d = double(dat(ch,:));
                tmp = sincFilter(d, ratio);
                validSamples = tmp(ntbuff+sampleRatio:sampleRatio:end);
                downData(ch,:) = validSamples;
            end
            dataSegments{nbChunks+1} = downData;
        else
            dataSegments{nbChunks+1} = [];
        end
    else
        dataSegments{nbChunks+1} = [];
    end

    fclose(fidI);

    % Concatenate all segments
    dataSegments = dataSegments(~cellfun('isempty',dataSegments));
    dataMatrix = cat(2, dataSegments{:});
    dataMatrix = double(dataMatrix); % ensure double

    % Construct timeVector and metaData
    numSamplesDown = size(dataMatrix,2);
    timeVector = (0:numSamplesDown-1)/outFs;

    metaData = struct();
    metaData.filePath = datFilePath;
    metaData.xmlPath = xmlFilePath;
    metaData.numChannels = nChannels;
    metaData.originalFs = inFs;
    metaData.samplingRate = outFs;
    metaData.numSamples = numSamplesDown;
    metaData.recordingTimeSec = numSamplesDown/outFs;
    metaData.lowPassFreq = lopass;

end
