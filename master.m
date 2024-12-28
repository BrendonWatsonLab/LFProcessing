% masterScript.m
%
% This master script demonstrates a basic workflow using the previously 
% defined functions. It assumes that you have enough RAM to load the entire
% dataset and that the .dat and .xml files are in the same directory with 
% the same base filename.
%
% Steps:
%   1) Load data and metadata using loadDatFile.
%   2) Plot raw voltage traces.
%   3) Generate spectrograms (now with functionality to select channels).
%   4) Plot PSDs.
%
% Author: Your Name
% Date: YYYY-MM-DD
% -------------------------------------------------------------

clear; close all; clc;

%% -------------------- File Setup --------------------
datFilePath = '/nfs/corexfs/MM-psych-watson-lab/EthanSCN/A471/A471.dat'; % Update with your actual path

%% -------------------- Load Data ---------------------
[dataMatrix, timeVector, metaData] = loadDatFile(datFilePath);

%% -------------------- Plot Raw Voltage Traces --------------------
% Plot raw traces for one channel. Default segments if none specified.
%channelNumber = 1;
%plotRawVoltageTraces(dataMatrix, timeVector, metaData, channelNumber);

%% -------------------- Generate Spectrograms --------------------
% Default: only channel 1, 24 hours or entire rec, 0-1000 Hz
%generateSpectrograms(dataMatrix, timeVector, metaData);

% If you want multiple channels, e.g. channels 1 and 2:
% generateSpectrograms(dataMatrix, timeVector, metaData, [1 2]);

% If you want all channels:
generateSpectrograms(dataMatrix, timeVector, metaData, 'all', [], [], [0 5]);

%% -------------------- Plot Channel PSDs --------------------
% Plot PSD for channel 1 in the default range [0 1000 Hz].
plotChannelPSDs(dataMatrix, metaData, 1, [0 5]);

% If you want multiple channels together:
% plotChannelPSDs(dataMatrix, metaData, [1 2 3 4]);

% End of masterScript.m
