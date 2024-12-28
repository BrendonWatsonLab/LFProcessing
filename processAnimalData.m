function processAnimalData(folderPath, basename, fs, freqBin, baselineRange, outputFile, varargin)
    % processAnimalData - Processes data for one animal, normalizes it to a baseline,
    %                     and saves a condensed summary in a .mat file.
    %
    % Syntax:
    %   processAnimalData(folderPath, basename, fs, freqBin, baselineRange, outputFile)
    %   processAnimalData(___, 'StartSample', startSample)
    %
    % Inputs:
    %   folderPath     - Path to the folder containing binary files for one animal.
    %   basename       - Base name of the binary files (e.g., 'A471').
    %   fs             - Original sampling frequency (Hz).
    %   freqBin        - Frequency bin (scalar or 2-element vector) for PSD averaging.
    %   baselineRange  - Range of time (in hours) to use as baseline [start, end].
    %   outputFile     - Path to save the condensed summary data.
    %   varargin       - Optional key-value pairs:
    %                       'StartSample': Starting sample for analysis (default = 1).
    %
    % Output:
    %   A .mat file containing normalized percent PSD values, mean, and SEM.

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

    % Process each channel to calculate normalized PSD
    allPSD = [];
    numChannels = length(binaryFiles);

    for i = 1:numChannels
        fileName = binaryFiles(i).name;
        filePath = fullfile(folderPath, fileName);
        fprintf('Processing file: %s\n', fileName);

        % Calculate PSD for this file
        [channelPSD, timeBins] = calculateAveragePSDFromBinary(filePath, fs, fs, freqBin, 1, 'StartSample', startSample);
        
        % Align time bins and normalize PSD to baseline
        baselineIdx = (timeBins >= baselineRange(1)) & (timeBins <= baselineRange(2));
        if ~any(baselineIdx)
            error('Baseline range does not overlap with available time bins.');
        end
        baselineMean = mean(channelPSD(baselineIdx));
        normalizedPSD = 100 * (channelPSD / baselineMean); % Normalize to baseline

        % Aggregate normalized PSD
        allPSD = [allPSD; normalizedPSD']; %#ok<AGROW>
    end

    % Calculate mean and SEM across channels
    meanPSD = mean(allPSD, 1);
    semPSD = std(allPSD, 0, 1) / sqrt(size(allPSD, 1));

    % Save condensed data
    save(outputFile, 'meanPSD', 'semPSD', 'timeBins', 'freqBin', 'baselineRange', '-v7.3');
    fprintf('Saved processed data to %s\n', outputFile);
end
