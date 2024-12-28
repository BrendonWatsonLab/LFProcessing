function generateSpectrogramsFromSplitFile(binFilePath, samplingRate, startTime_s, duration_s, freqRange)
    %% Defaults
    if ~exist('startTime_s','var') || isempty(startTime_s)
        startTime_s = 0;
    end
    if ~exist('duration_s','var') || isempty(duration_s)
        duration_s = inf; % Default to entire file
    end
    if ~exist('freqRange','var') || isempty(freqRange)
        freqRange = [0 100];
    end

    %% File Information
    fileID = fopen(binFilePath, 'r');
    fseek(fileID, 0, 'eof');
    totalSamples = ftell(fileID) / 8; % Assuming double precision (8 bytes per sample)
    fclose(fileID);

    totalDuration = totalSamples / samplingRate;

    % Adjust duration if it exceeds the file length
    if duration_s == inf || (startTime_s + duration_s) > totalDuration
        duration_s = totalDuration - startTime_s;
    end

    %% Downsample if Necessary
    maxFreq = freqRange(end);
    targetFs = max(2 * maxFreq, 2 * (maxFreq + 5)); % Nyquist frequency rule
    downsampleFactor = ceil(samplingRate / targetFs);

    if downsampleFactor > 1
        fprintf('Downsampling by a factor of %d to reduce sampling rate to %.2f Hz.\n', downsampleFactor, samplingRate / downsampleFactor);
    end

    %% Adjusted Chunk Parameters
    minSamplesForWindow = 2^10; % Minimum samples to support at least two windows
    chunkSize = max(round((samplingRate / downsampleFactor) * 2), minSamplesForWindow);
    numChunks = ceil((duration_s * samplingRate / downsampleFactor) / chunkSize);

    %% Spectrogram Parameters
    winLength = max(round((samplingRate / downsampleFactor) * 0.5), 128); % 0.5-second window
    overlap = round(winLength / 2);
    w = 0.5 - 0.5 * cos(2 * pi * (0:winLength-1)' / (winLength-1)); % Hann window manually implemented
    windowPower = sum(w.^2);
    nfft = max(2^nextpow2(winLength), 256); % Ensure adequate resolution

    %% Initialize Spectrogram Accumulators
    S_accum = [];
    T_accum = [];
    F_accum = [];

    %% Process in Chunks
    for chunkIdx = 1:numChunks
        chunkStartSample = (chunkIdx - 1) * chunkSize * downsampleFactor + max(1, floor(startTime_s * samplingRate) + 1);
        chunkEndSample = min(chunkStartSample + chunkSize * downsampleFactor - 1, totalSamples);

        % Read chunk
        fileID = fopen(binFilePath, 'r');
        fseek(fileID, (chunkStartSample - 1) * 8, 'bof');
        chunkData = fread(fileID, [1, chunkEndSample - chunkStartSample + 1], 'double');
        fclose(fileID);

        % Downsample the chunk if necessary
        if downsampleFactor > 1
            chunkData = downsample(chunkData, downsampleFactor);
        end

        % Skip chunks that are too small
        if length(chunkData) < winLength
            warning('Chunk too small for spectrogram computation. Skipping chunk %d.', chunkIdx);
            continue;
        end

        % Compute Spectrogram for Chunk
        [S, F, T] = spectrogram(chunkData, w, overlap, nfft, samplingRate / downsampleFactor);
        P = abs(S).^2 / (windowPower * (samplingRate / downsampleFactor)); % PSD in V^2/Hz

        % Subset frequencies
        validFreqIdx = (F >= freqRange(1) & F <= freqRange(2));
        if ~any(validFreqIdx)
            error('No frequency bins found in the specified range [%.2f %.2f] Hz. Increase nfft or adjust freqRange.', freqRange(1), freqRange(2));
        end

        Fsub = F(validFreqIdx);
        Psub = P(validFreqIdx, :);

        % Accumulate Results
        S_accum = [S_accum, Psub];
        T_accum = [T_accum, T + (chunkIdx - 1) * (chunkSize / (samplingRate / downsampleFactor))];
        F_accum = Fsub;
    end

    %% Plot and Save Spectrogram
    fig = figure('Visible', 'off', 'Color', 'white');
    imagesc(T_accum, F_accum, 10*log10(S_accum));
    axis xy;
    colormap(parula);
    colorbar;
    xlabel('Time (s)', 'FontName', 'Arial', 'FontSize', 14);
    ylabel('Frequency (Hz)', 'FontName', 'Arial', 'FontSize', 14);
    title(sprintf('Spectrogram: %.1f-%.1f Hz', freqRange(1), freqRange(2)), 'FontName', 'Arial', 'FontSize', 16);
    fig.Position = [100, 100, 1200, 200];

    % Save the figure
    [~, baseName, ~] = fileparts(binFilePath);
    pngFileName = sprintf('%s_spectrogram_%.0fTo%.0fHz.png', baseName, freqRange(1), freqRange(2));
    exportgraphics(fig, pngFileName, 'Resolution', 300);
    close(fig);

    fprintf('Spectrogram saved to %s\n', pngFileName);
end
