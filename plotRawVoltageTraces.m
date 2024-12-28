function plotRawVoltageTraces(dataMatrix, timeVector, metaData, channelNumber, ...
                              dailyStartTime_s, hourlyStartTime_s, ...
                              minutelyStartTime_s, secondlyStartTime_s)
% plotRawVoltageTraces - Generate and save PNGs of raw voltage traces from a specified channel.
%
%   This function creates PNG images of voltage vs. time for a specified channel
%   over four different durations: 
%       1 day (86400 s), 
%       1 hour (3600 s),
%       1 minute (60 s), 
%       1 second (1 s).
%
%   It also automatically adjusts the vertical scale (y-axis limits) to fit
%   the data for each plotted segment, preventing the trace from appearing 
%   as a flattened line.
%
%   If start times are not provided, the function will select windows from 
%   the middle of the recording (or as close as possible if shorter).
%
%   The resulting figures are saved as PNG files. The figures are formatted 
%   as long horizontal rectangles, with a basic y-axis showing voltage and 
%   no visible x-axis. The voltage trace is plotted in dark green.
%
%   Inputs:
%       dataMatrix      - [nChannels x nSamples] matrix of raw data.
%       timeVector      - [1 x nSamples] time vector (in seconds).
%       metaData        - Struct with fields:
%                           * samplingRate
%                           * numChannels
%                           * recordingTimeSec
%       channelNumber   - (integer) index of the channel to plot.
%
%       dailyStartTime_s    - (optional) start time for the 1-day plot.
%       hourlyStartTime_s   - (optional) start time for the 1-hour plot.
%       minutelyStartTime_s - (optional) start time for the 1-minute plot.
%       secondlyStartTime_s - (optional) start time for the 1-second plot.
%
%   Outputs:
%       No direct outputs. The function saves PNG figures.
%
%   Example:
%       plotRawVoltageTraces(dataMatrix, timeVector, metaData, 1);
%
%   Author: Your Name
%   Date: YYYY-MM-DD
%   -------------------------------------------------------------

    %% -------------------- Parameter & Default Handling --------------------
    % Durations (in seconds)
    dayDur    = 86400; % 24 hours
    hourDur   = 3600;  % 1 hour
    minuteDur = 60;    % 1 minute
    secondDur = 1;     % 1 second

    totalDuration = metaData.recordingTimeSec;

    % If start times are not provided, pick from the middle if possible
    if ~exist('dailyStartTime_s','var') || isempty(dailyStartTime_s)
        dailyStartTime_s = max(0, totalDuration/2 - dayDur/2);
    end
    if ~exist('hourlyStartTime_s','var') || isempty(hourlyStartTime_s)
        hourlyStartTime_s = max(0, totalDuration/2 - hourDur/2);
    end
    if ~exist('minutelyStartTime_s','var') || isempty(minutelyStartTime_s)
        minutelyStartTime_s = max(0, totalDuration/2 - minuteDur/2);
    end
    if ~exist('secondlyStartTime_s','var') || isempty(secondlyStartTime_s)
        secondlyStartTime_s = max(0, totalDuration/2 - secondDur/2);
    end

    %% -------------------- Validate Channel Number --------------------
    if channelNumber < 1 || channelNumber > metaData.numChannels
        error('Invalid channelNumber. Must be between 1 and %d.', metaData.numChannels);
    end

    %% -------------------- Extract Channel Data --------------------
    channelData = dataMatrix(channelNumber, :);

    %% -------------------- Helper Function for Plotting --------------------
    function plotAndSavePNG(startTime, duration, filenameSuffix)
        endTime = startTime + duration;
        if endTime > totalDuration
            endTime = totalDuration;
        end

        fs = metaData.samplingRate;
        startIndex = max(1, floor(startTime * fs) + 1);
        endIndex   = min(length(timeVector), ceil(endTime * fs));
        
        tChunk = timeVector(startIndex:endIndex);
        vChunk = channelData(startIndex:endIndex);

        % Create figure
        fig = figure('Visible','off','Color','white');

        % Plot with dark green color
        plot(tChunk, vChunk, 'Color', [0 0.5 0], 'LineWidth', 1);

        % Format axes
        ax = gca;
        ax.FontName = 'Arial';
        ax.FontSize = 12;
        ax.Box = 'off'; 
        ax.XColor = 'none'; % no visible x-axis
        ax.XTick = [];
        ylabel('Voltage (µV)');

        % Adjust y-limits to fit data
        dataMin = min(vChunk);
        dataMax = max(vChunk);
        % Add a small margin:
        margin = 0.05 * (dataMax - dataMin);
        if margin == 0
            margin = 1; % in case data is flat
        end
        ylim([dataMin - margin, dataMax + margin]);

        % Adjust figure size (Width x Height in pixels)
        fig.Position = [100 100 1200 200];

        % Save as PNG
        pngFileName = sprintf('channel%d_%s.png', channelNumber, filenameSuffix);
        exportgraphics(fig, pngFileName,'Resolution',300);

        % Close the figure
        close(fig);
    end

    %% -------------------- Plot & Save Each Duration --------------------
    % If the requested duration is longer than the recording, plot what is available.

    % 1 Day
    if dayDur <= totalDuration
        plotAndSavePNG(dailyStartTime_s, dayDur, '1day');
    else
        plotAndSavePNG(0, totalDuration, 'entireRecordingFor1dayRequest');
    end

    % 1 Hour
    if hourDur <= totalDuration
        plotAndSavePNG(hourlyStartTime_s, hourDur, '1hour');
    else
        plotAndSavePNG(0, totalDuration, 'entireRecordingFor1hourRequest');
    end

    % 1 Minute
    if minuteDur <= totalDuration
        plotAndSavePNG(minutelyStartTime_s, minuteDur, '1minute');
    else
        plotAndSavePNG(0, totalDuration, 'entireRecordingFor1minuteRequest');
    end

    % 1 Second
    if secondDur <= totalDuration
        plotAndSavePNG(secondlyStartTime_s, secondDur, '1second');
    else
        plotAndSavePNG(0, totalDuration, 'entireRecordingFor1secondRequest');
    end

end

