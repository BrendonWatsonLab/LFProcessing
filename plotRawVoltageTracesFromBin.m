function plotRawVoltageTracesFromBin(binFilePath, xmlFilePath, dailyStartTime_s, hourlyStartTime_s, minutelyStartTime_s, secondlyStartTime_s)
% plotRawVoltageTracesFromBin - Generate and save PNGs of raw voltage traces from a .bin file.
%
%   This function creates PNG images of voltage vs. time for a specified channel
%   over four different durations:
%       1 day (86400 s),
%       1 hour (3600 s),
%       1 minute (60 s),
%       1 second (1 s).
%
%   It reads data from a binary file for the specified channel and obtains metadata
%   (e.g., sampling rate) from an accompanying XML file. It automatically adjusts
%   the vertical scale (y-axis limits) to fit the data for each plotted segment, and
%   saves the resulting figures as PNG files.
%
%   Inputs:
%       binFilePath        - String, path to the .bin file for the specified channel.
%       xmlFilePath        - String, path to the .xml file containing metadata.
%       dailyStartTime_s   - (optional) start time for the 1-day plot.
%       hourlyStartTime_s  - (optional) start time for the 1-hour plot.
%       minutelyStartTime_s- (optional) start time for the 1-minute plot.
%       secondlyStartTime_s- (optional) start time for the 1-second plot.
%
%   Outputs:
%       No direct outputs. The function saves PNG figures.
%
%   Example:
%       plotRawVoltageTracesFromBin('channel1.bin', 'metadata.xml');
%
%   Author: Your Name
%   Date: YYYY-MM-DD
%   -------------------------------------------------------------

    %% -------------------- Parse XML Metadata --------------------
    % Read XML file and extract metadata
    xmlDoc = xmlread(xmlFilePath);
    samplingRateNode = xmlDoc.getElementsByTagName('samplingRate').item(0);
    samplingRate = str2double(samplingRateNode.getTextContent());

    % Open the binary file and read all data
    fileID = fopen(binFilePath, 'r');
    if fileID == -1
        error('Failed to open file: %s', binFilePath);
    end

    channelData = fread(fileID, 'float32'); % Assuming float32 data format
    fclose(fileID);

    % Determine recording duration from the length of channelData
    totalDuration = length(channelData) / samplingRate;

    %% -------------------- Parameter & Default Handling --------------------
    % Durations (in seconds)
    dayDur    = 86400; % 24 hours
    hourDur   = 3600;  % 1 hour
    minuteDur = 60;    % 1 minute
    secondDur = 1;     % 1 second

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

    %% -------------------- Helper Function for Plotting --------------------
    function plotAndSavePNG(startTime, duration, filenameSuffix)
        endTime = startTime + duration;
        if endTime > totalDuration
            endTime = totalDuration;
        end

        fs = samplingRate;
        startIndex = max(1, floor(startTime * fs) + 1);
        endIndex   = min(length(channelData), ceil(endTime * fs));

        tChunk = (startIndex:endIndex) / fs; % Convert sample indices to time
        vChunk = channelData(startIndex:endIndex);

        % Create figure
        fig = figure('Visible','off','Color','white');

        % Plot with dark green color
        plot(tChunk, vChunk, 'Color', [0 0.5 0], 'LineWidth', 0.1);

        % Format axes
        ax = gca;
        ax.FontName = 'Arial';
        ax.FontSize = 12;
        ax.Box = 'off'; 
        ax.XColor = 'none'; % no visible x-axis
        ax.XTick = [];
        ylabel('Voltage (\muV)');

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
        pngFileName = sprintf('channel_%s.png', filenameSuffix);
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
