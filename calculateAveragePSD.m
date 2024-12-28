function calculatePSDForAllChannels(folderPath, basename, fs, freqBin, varargin)
    % calculatePSDForAllChannels - Computes average PSD and SEM for all binary files
    %                               matching the basename and generates PNG plots.
    %
    % Syntax:
    %   calculatePSDForAllChannels(folderPath, basename, fs, freqBin)
    %   calculatePSDForAllChannels(folderPath, basename, fs, freqBin, 'StartSample', startSample)
    %
    % Inputs:
    %   folderPath - Path to the folder containing the binary files.
    %   basename   - Base name of the binary files (e.g., 'A471').
    %   fs         - Original sampling frequency (Hz).
    %   freqBin    - Frequency bin (scalar or 2-element vector) for PSD averaging.
    %   varargin   - Optional key-value pairs:
    %                   'StartSample': Starting sample for analysis (default = 1).
    %
    % Outputs:
    %   Generates PNG plots saved in the same folder.

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'StartSample', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    startSample = p.Results.StartSample;

    % Get list of binary files matching the pattern
    filePattern = fullfile(folderPath, sprintf('%s_channel*.bin', basename));
    binaryFiles = dir(filePattern);
    
    if isempty(binaryFiles)
        error('No binary files found matching pattern: %s', filePattern);
    end
    
    % Frequency bin adjustment
    if numel(freqBin) == 1
        freqBin = [freqBin freqBin]; % Treat scalar as a single bin
    elseif numel(freqBin) ~= 2
        error('Frequency bin must be a scalar or a 2-element vector.');
    end
    
    % Determine target sampling frequency for downsampling
    maxFreq = max(freqBin);
    targetFs = min(fs, 2 * maxFreq); % Ensure at least Nyquist frequency
    downsampleFactor = round(fs / targetFs);

    % Process each file
    for i = 1:length(binaryFiles)
        fileName = binaryFiles(i).name;
        filePath = fullfile(folderPath, fileName);
        
        % Extract channel number from file name
        channelMatch = regexp(fileName, sprintf('%s_channel(\\d+)', basename), 'tokens');
        if isempty(channelMatch)
            warning('Skipping file %s: does not match expected format.', fileName);
            continue;
        end
        channelNum = str2double(channelMatch{1}{1});
        
        % Process the binary file
        fprintf('Processing file: %s (Channel %d)\n', fileName, channelNum);
        [avgPSD, semPSD] = calculateAveragePSDFromBinary(filePath, fs, targetFs, freqBin, downsampleFactor, 'StartSample', startSample);

        % Generate plot with SEM error bars
        timeBins = (1:length(avgPSD)) * (3600 / 3600); % Time in hours
        figure;
        errorbar(timeBins, avgPSD, semPSD, '-o', 'LineWidth', 2);
        xlabel('Time (hours)');
        ylabel('Average PSD (Power/Frequency)');
        title(sprintf('Channel %d: Average PSD with SEM for %d-%d Hz', channelNum, freqBin(1), freqBin(2)));
        grid on;

        % Save plot as PNG
        outputFileName = sprintf('%s_channel%d_PSD.png', basename, channelNum);
        outputFilePath = fullfile(folderPath, outputFileName);
        saveas(gcf, outputFilePath);
        fprintf('Saved plot to %s\n', outputFilePath);

        % Close figure to free memory
        close;
    end
end

function [avgPSD, semPSD] = calculateAveragePSDFromBinary(filePath, originalFs, targetFs, freqBin, downsampleFactor, varargin)
    % Sub-function to calculate PSD and SEM for a single binary file with downsampling

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'StartSample', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    startSample = p.Results.StartSample;

    % Constants
    oneHourSamples = targetFs * 3600; % Samples in one hour after downsampling
    window = round(targetFs * 2); % 2-second window for Welch
    overlap = round(window / 2); % 50% overlap
    nfft = max(256, 2^nextpow2(window)); % NFFT for FFT
    
    % Open the binary file
    fid = fopen(filePath, 'r');
    if fid == -1
        error('Failed to open the file: %s', filePath);
    end
    
    % Get file size and calculate total samples
    fseek(fid, 0, 'eof');
    fileSize = ftell(fid);
    totalSamples = fileSize / 8; % Assuming double-precision (8 bytes/sample)
    fclose(fid);
    
    % Adjust start sample
    startSample = max(1, min(startSample, totalSamples));

    % Initialize the output
    avgPSD = [];
    semPSD = [];

    % Process in 1-hour bins
    for binStart = startSample:oneHourSamples * downsampleFactor:totalSamples
        binEnd = min(binStart + oneHourSamples * downsampleFactor - 1, totalSamples);
        if binEnd <= binStart
            break;
        end
        
        % Read the data for the current bin
        fid = fopen(filePath, 'r');
        fseek(fid, (binStart - 1) * 8, 'bof'); % Skip to bin start
        segment = fread(fid, binEnd - binStart + 1, 'double');
        fclose(fid);
        
        % Downsample the data
        segment = decimate(segment, downsampleFactor);
        
        % Calculate PSD using Welch's method
        [Pxx, f] = pwelch(segment, window, overlap, nfft, targetFs);

        % Identify the frequency indices within the specified bin
        freqIdx = (f >= freqBin(1)) & (f <= freqBin(2));
        if ~any(freqIdx)
            error('Specified frequency bin is out of range for the data.');
        end

        % Compute average PSD and SEM in the frequency bin
        psdValues = Pxx(freqIdx);
        avgPSD = [avgPSD; mean(psdValues)]; %#ok<AGROW>
        semPSD = [semPSD; std(psdValues) / sqrt(length(psdValues))]; %#ok<AGROW>
    end
end
