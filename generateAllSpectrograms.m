function generateAllSpectrograms(baseName, samplingRate, startTime_s, duration_s, freqRange)
% generateAllSpectrograms - Generate spectrograms for all channel files.
% 
%   This function iterates over all channel files with names formatted as
%   "<baseName>_channelX.bin" where X is the channel number (1 to 30),
%   and generates spectrograms for each channel using the
%   `generateSpectrogramsFromSplitFile` function.
% 
%   Inputs:
%       baseName      - Base name of the channel files (e.g., 'A471').
%       samplingRate  - Sampling rate of the data (Hz).
%       startTime_s   - Start time in seconds. Default: 0.
%       duration_s    - Duration in seconds. Default: entire file.
%       freqRange     - [fStart fEnd] Hz. Default: [0 100].
% 
%   Outputs:
%       Generates and saves spectrogram images for each channel file.
% 
%   Example:
%       generateAllSpectrograms('A471', 20000, 0, 60, [0 5]);
% 
%   Author: Your Name
%   Date: YYYY-MM-DD
%   -------------------------------------------------------------

    %% Defaults
    if ~exist('startTime_s', 'var') || isempty(startTime_s)
        startTime_s = 0;
    end
    if ~exist('duration_s', 'var') || isempty(duration_s)
        duration_s = inf; % Default to entire file
    end
    if ~exist('freqRange', 'var') || isempty(freqRange)
        freqRange = [0 100];
    end

    %% Iterate over all channels (1 to 30)
    for channelNum = 1:30
        % Construct the channel file name
        binFilePath = sprintf('%s_channel%d.bin', baseName, channelNum);

        % Check if the file exists
        if ~isfile(binFilePath)
            warning('File %s does not exist. Skipping...', binFilePath);
            continue;
        end

        % Call the spectrogram generation function
        fprintf('Processing channel %d...', channelNum);
        try
            generateSpectrogramsFromSplitFile(binFilePath, samplingRate, startTime_s, duration_s, freqRange);
        catch ME
            fprintf('Error processing channel %d: %s', channelNum, ME.message);
        end
    end

    fprintf('Spectrogram generation completed for all channels.');
end
