function generateSpectrograms(dataMatrix, timeVector, metaData, channels, startTime_s, duration_s, freqRange, downsampleFactor)
% generateSpectrograms - Generate and save spectrogram images with filtering, optional downsampling,
%                        and proper frequency subsetting.
%
%   This version:
%   - Uses a large nfft to ensure fine frequency resolution in the low-frequency range.
%   - Subsets frequencies after computing the spectrogram.
%   - If no frequencies are found in the requested range, it raises an error.
%
%   Inputs:
%       dataMatrix      - [nChannels x nSamples] raw data (volts)
%       timeVector      - [1 x nSamples] time vector (seconds)
%       metaData        - struct with:
%                           * samplingRate
%                           * numChannels
%                           * recordingTimeSec
%       channels        - 'all', single channel number, or array of channels. Default: 1
%       startTime_s     - start time in seconds. Default: 0
%       duration_s      - duration in seconds. Default: 86400 (24 hours)
%       freqRange       - [fStart fEnd] Hz. Default: [0 100]
%       downsampleFactor- integer factor for downsampling. Default: 1 (no downsampling)
%
%   Example:
%       generateSpectrograms(dataMatrix, timeVector, metaData, 'all', [], [], [0 5]);
%
%   Author: Your Name
%   Date: YYYY-MM-DD
%   -------------------------------------------------------------

    %% Defaults
    if ~exist('channels','var') || isempty(channels)
        channels = 1; 
    end
    if ischar(channels) && strcmpi(channels, 'all')
        channels = 1:metaData.numChannels;
    end
    if ~exist('startTime_s','var') || isempty(startTime_s)
        startTime_s = 0;
    end
    if ~exist('duration_s','var') || isempty(duration_s)
        duration_s = 86400; % 24h
    end
    if ~exist('freqRange','var') || isempty(freqRange)
        freqRange = [0 100];
    end
    if ~exist('downsampleFactor','var') || isempty(downsampleFactor)
        downsampleFactor = 1;
    end

    fs = metaData.samplingRate;
    totalDuration = metaData.recordingTimeSec;

    % Adjust duration if needed
    if startTime_s < 0
        startTime_s = 0;
    end
    if startTime_s + duration_s > totalDuration
        duration_s = totalDuration - startTime_s;
    end

    startIdx = max(1, floor(startTime_s * fs) + 1);
    endIdx   = min(length(timeVector), startIdx + floor(duration_s * fs) - 1);

    segmentData = dataMatrix(:, startIdx:endIdx);
    segmentTime = timeVector(startIdx:endIdx);

    %% Downsampling
    if downsampleFactor > 1
        newFs = fs / downsampleFactor;
        if newFs < 2 * freqRange(end)
            warning('Downsampling too aggressive, Nyquist < max freq. Consider adjusting.');
        end
        segmentData = downsample(segmentData.', downsampleFactor).';
        segmentTime = downsample(segmentTime, downsampleFactor);
        fs = newFs;
    end

    fStart = freqRange(1);
    fEnd   = freqRange(end);
    if fEnd > fs/2
        warning('Upper frequency > Nyquist. Adjusting to Nyquist.');
        fEnd = fs/2;
    end

    %% Filter Design
    filterOrder = 4;
    if fStart <= 1e-10
        % If starting at 0 Hz, design a lowpass
        Wn = fEnd/(fs/2);
        if Wn > 1, Wn = 1; end
        [b,a] = butter(filterOrder, Wn, 'low');
    else
        % Bandpass
        Wn = [fStart/(fs/2), fEnd/(fs/2)];
        if Wn(1) >= Wn(2)
            error('Invalid band range.');
        end
        [b,a] = butter(filterOrder, Wn, 'bandpass');
    end

    %% Spectrogram Parameters
    % Increase nfft significantly for fine frequency resolution
    % For 0-5 Hz at 20 kHz, use large nfft:
    nfft = 131072; % Frequency spacing ~ fs/nfft ~ 20000/131072 ~ 0.15 Hz
    winLength = round(fs * 1.0); % 1 second
    overlap   = round(winLength/2);
    w = bartlett(winLength);
    windowPower = sum(w.^2);

    %% Gaussian Smoothing Kernel
    sigma = 1;
    kernelSize = 7;
    x = linspace(-3,3,kernelSize);
    gaussKernel = exp(-x.^2/(2*sigma^2));
    gaussKernel = gaussKernel / sum(gaussKernel);

    %% Compute & Save Spectrograms
    for ch = channels
        chData = filtfilt(b, a, segmentData(ch,:));

        [S,F,T] = spectrogram(chData, w, overlap, nfft, fs);
        P = abs(S).^2 / (windowPower * fs); % PSD in V^2/Hz

        % Smooth along frequency
        P_smoothed = conv2(P, gaussKernel', 'same');

        % Subset frequencies
        validFreqIdx = (F >= fStart & F <= fEnd);
        if ~any(validFreqIdx)
            error('No frequency bins found in the specified range [%.2f %.2f] Hz. Increase nfft or adjust freqRange.', fStart, fEnd);
        end

        Fsub = F(validFreqIdx);
        Psub = P_smoothed(validFreqIdx,:);

        fig = figure('Visible','off','Color','white');
        imagesc(T, Fsub, 10*log10(Psub));
        axis xy;
        colormap(parula);
        c = colorbar;
        c.Label.String = 'Power (dB V^2/Hz)';
        c.FontName = 'Arial';
        c.FontSize = 12;

        xlabel('Time (s)', 'FontName', 'Arial', 'FontSize', 14);
        ylabel('Frequency (Hz)', 'FontName', 'Arial', 'FontSize', 14);

        ax = gca;
        ax.FontName = 'Arial';
        ax.FontSize = 12;
        ax.Box = 'off';

        title(sprintf('Ch %d: %.1f-%.1f Hz (Downx%d, %ds to %ds)', ch, fStart, fEnd, downsampleFactor, startTime_s, startTime_s+duration_s), ...
            'FontName', 'Arial', 'FontSize', 16);

        fig.Position = [100 100 1200 600];

        pngFileName = sprintf('channel%d_spectrogram_%.0fsecTo%.0fsec_%.1fTo%.1fHz_downx%d.png', ...
                              ch, startTime_s, startTime_s+duration_s, fStart, fEnd, downsampleFactor);
        exportgraphics(fig, pngFileName,'Resolution',300);
        close(fig);
    end

end