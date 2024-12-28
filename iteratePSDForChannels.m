function iteratePSDForChannels(folderPath, baseName, fs, freqRange)
    % iteratePSDForChannels: Function to compute PSD for each channel file in a folder
    %
    % Parameters:
    %   folderPath - Path to the folder containing the .bin files (string)
    %   baseName - Base name of the files (string, e.g., 'A471')
    %   fs - Sampling frequency in Hz (scalar)
    %   freqRange - Frequency range to plot [minFreq, maxFreq] (1x2 vector)

    % Get a list of files matching the base name and channel pattern
    filePattern = fullfile(folderPath, [baseName, '_channel*.bin']);
    channelFiles = dir(filePattern);

    if isempty(channelFiles)
        error('No files found matching the pattern %s', filePattern);
    end

    % Iterate through each file and compute PSD
    for i = 1:length(channelFiles)
        channelFilePath = fullfile(folderPath, channelFiles(i).name);
        fprintf('Processing file: %s\n', channelFiles(i).name);

        % Call plotPSDfromBin for each file
        try
            plotPSDfromBin(channelFilePath, fs, freqRange);
        catch ME
            fprintf('Error processing file %s: %s\n', channelFiles(i).name, ME.message);
        end
    end

    fprintf('Processed %d files in folder %s\n', length(channelFiles), folderPath);
end
