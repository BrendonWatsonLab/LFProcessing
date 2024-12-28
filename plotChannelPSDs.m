function plotChannelPSDs(dataMatrix, metaData, channelNumbers, freqRange)
% plotChannelPSDs - Plot PSD (Power Spectral Density) for one or multiple channels.
%
%   This function computes and plots the PSD for the specified channel(s) 
%   from the dataMatrix. Frequency is on the x-axis, and PSD on the y-axis.
%
%   If multiple channels are provided, it plots all PSDs on the same graph, 
%   with distinct line colors for each channel.
%
%   Inputs:
%       dataMatrix    - [nChannels x nSamples] raw data matrix
%       metaData      - Struct with fields:
%                         * samplingRate: sampling frequency (Hz)
%                         * numChannels
%                         * recordingTimeSec (optional)
%       channelNumbers - Integer or array of integers specifying which channels 
%                        to plot. (1-based indexing)
%       freqRange     - [fStart fEnd] frequency range in Hz for plotting PSD.
%                       If not provided, defaults to [0 1000].
%
%   Outputs:
%       No direct outputs. This function creates and saves a PNG figure.
%
%   Notes:
%       * The PSD is estimated using Welch's method (pwelch).
%       * For multiple channels, each channel's line is plotted in a different color.
%       * PSD units will be in power per Hz (dB/Hz or linear scale depending on plotting).
%         Here we choose to plot in decibels to be consistent.
%
%   Example:
%       plotChannelPSDs(dataMatrix, metaData, 1);
%       plotChannelPSDs(dataMatrix, metaData, [1 2 3], [0 500]);
%
%   Author: Your Name
%   Date: YYYY-MM-DD
%   -------------------------------------------------------------

    %% -------------------- Default Parameters --------------------
    if ~exist('freqRange','var') || isempty(freqRange)
        freqRange = [0 1000]; % default frequency range
    end

    fs = metaData.samplingRate;

    % Validate channel numbers
    if any(channelNumbers < 1) || any(channelNumbers > metaData.numChannels)
        error('Invalid channelNumbers. Must be between 1 and %d.', metaData.numChannels);
    end

    %% -------------------- Compute PSD for Each Channel --------------------
    % Use pwelch to compute PSD. Adjust parameters as needed.
    % We'll use a standard window (e.g., 1024 samples), 50% overlap, and 1024-point FFT.
    winLength = 1024;
    overlap   = 512;
    nfft      = 1024;

    allPSD = {};   % cell to store PSDs
    allFreqs = {}; % cell to store frequency vectors

    for ch = channelNumbers
        [Pxx, F] = pwelch(dataMatrix(ch,:), winLength, overlap, nfft, fs);
        % Convert to dB
        Pxx_dB = 10*log10(Pxx);
        
        allPSD{end+1} = Pxx_dB;
        allFreqs{end+1} = F;
    end

    %% -------------------- Select Frequency Range --------------------
    % We'll trim the frequencies to the requested freqRange.
    fMin = freqRange(1);
    fMax = freqRange(end);

    % Assume allFreqs{1} covers the full range since all PSD computations are identical in method.
    freqIdx = (allFreqs{1} >= fMin & allFreqs{1} <= fMax);
    plotF = allFreqs{1}(freqIdx);

    %% -------------------- Plotting --------------------
    fig = figure('Visible','off','Color','white');
    ax = gca;
    hold on;

    % Set up a colormap for multiple channels
    numCh = length(channelNumbers);
    cmap = lines(numCh);

    for i = 1:numCh
        plotP = allPSD{i}(freqIdx);
        plot(plotF, plotP, 'LineWidth', 1.5, 'Color', cmap(i,:));
    end

    % Formatting
    xlabel('Frequency (Hz)', 'FontName', 'Arial', 'FontSize', 14);
    ylabel('PSD (dB/Hz)', 'FontName', 'Arial', 'FontSize', 14);
    titleStr = 'PSD for ';
    if numCh == 1
        titleStr = [titleStr sprintf('Channel %d', channelNumbers)];
    else
        titleStr = [titleStr 'Channels ' num2str(channelNumbers)];
    end
    title(titleStr, 'FontName', 'Arial', 'FontSize', 16);

    ax.FontName = 'Arial';
    ax.FontSize = 12;
    ax.Box = 'off';
    grid on;

    if numCh > 1
        legendStrings = arrayfun(@(c) sprintf('Ch %d', c), channelNumbers, 'UniformOutput', false);
        legend(legendStrings, 'Location','best', 'Box','off');
    end

    % Adjust figure size (width x height in pixels)
    fig.Position = [100 100 1200 600];

    %% -------------------- Save the Figure --------------------
    % Construct filename
    if numCh == 1
        pngFileName = sprintf('channel%d_psd_%.0fTo%.0fHz.png', channelNumbers, fMin, fMax);
    else
        chStr = strjoin(arrayfun(@(c) sprintf('%d', c), channelNumbers, 'UniformOutput', false), '_');
        pngFileName = sprintf('channels_%s_psd_%.0fTo%.0fHz.png', chStr, fMin, fMax);
    end

    exportgraphics(fig, pngFileName, 'Resolution',300);

    % Close figure
    close(fig);

end
