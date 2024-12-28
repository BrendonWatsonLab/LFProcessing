function plotPSDfromBin(binFilePath, fs, freqRange)
    % plotPSDfromBin: Function to plot and save average PSD vs frequency from a single-channel .bin file
    %
    % Parameters:
    %   binFilePath - Path to the .bin file (string)
    %   fs - Sampling frequency in Hz (scalar)
    %   freqRange - Frequency range to plot [minFreq, maxFreq] (1x2 vector)

    % Read .bin file (assuming data is stored as double precision floats)
    fileID = fopen(binFilePath, 'rb');
    rawData = fread(fileID, 'double');
    fclose(fileID);

    % Check if rawData is empty
    if isempty(rawData)
        error('The .bin file appears to be empty or improperly read.');
    end

    % Auto downscale to match the specified frequency range
    minFreq = freqRange(1);
    maxFreq = freqRange(2);
    targetFs = max(2 * maxFreq, 2 * (maxFreq + 5)); % Nyquist frequency rule

    downsampleFactor = max(1, ceil(fs / targetFs));
    if downsampleFactor > 1
        rawData = downsample(rawData, downsampleFactor);
        fs = fs / downsampleFactor;
    end

    % Ensure rawData is a vector (single-channel assumption)
    numSamples = length(rawData);

    % Set segment parameters similar to spectrogram
    winLength = max(round((fs) * 0.5), 128); % 0.5-second window
    overlap = round(winLength / 2);
    w = 0.5 - 0.5 * cos(2 * pi * (0:winLength-1)' / (winLength-1)); % Hann window
    windowPower = sum(w.^2);
    nfft = max(2^nextpow2(winLength), 256);

    % Compute PSD using Welch's method with spectrogram-like parameters
    [pxx, f] = pwelch(rawData, w, overlap, nfft, fs);
    pxx = pxx / windowPower; % Normalize PSD similar to spectrogram

    % Limit frequency range
    freqMask = (f >= freqRange(1)) & (f <= freqRange(2));
    f = f(freqMask);
    pxx = pxx(freqMask);

    % Check if f and pxx are empty after masking
    if isempty(f) || isempty(pxx)
        error('The specified frequency range does not overlap with the computed spectrum.');
    end

    % Plot average PSD
    figure;
    plot(f, 10 * log10(pxx), 'b');
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density (dB/Hz)');
    title('Average Power Spectral Density vs Frequency');
    grid on;

    % Save the plot as a PNG file
    [~, fileName, ~] = fileparts(binFilePath);
    saveFileName = sprintf('%s_PSD.png', fileName);
    saveas(gcf, saveFileName);
    fprintf('Saved PSD plot to %s\n', saveFileName);

    % Close the figure to prevent clutter
    close(gcf);
end
