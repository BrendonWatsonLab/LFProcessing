function combineAnimalSummaries(summaryFiles, outputFile)
    % combineAnimalSummaries - Combines multiple summary files to calculate the
    %                          final average PSD and SEM across all animals.
    %
    % Syntax:
    %   combineAnimalSummaries(summaryFiles, outputFile)
    %
    % Inputs:
    %   summaryFiles - Cell array of paths to .mat files containing animal summaries.
    %   outputFile   - Path to save the final plot and combined data.
    %
    % Output:
    %   A PNG plot and a .mat file with combined mean PSD and SEM.

    % Initialize storage for aggregated data
    aggregatedMeanPSD = [];
    aggregatedSEMPSD = [];
    combinedTimeBins = [];
    
    for i = 1:length(summaryFiles)
        % Load each summary file
        data = load(summaryFiles{i});
        if isempty(aggregatedMeanPSD)
            combinedTimeBins = data.timeBins; % Use time bins from the first file
        else
            % Ensure time bins match across files
            if any(data.timeBins ~= combinedTimeBins)
                error('Time bins do not align across summary files.');
            end
        end
        aggregatedMeanPSD = [aggregatedMeanPSD; data.meanPSD]; %#ok<AGROW>
        aggregatedSEMPSD = [aggregatedSEMPSD; data.semPSD]; %#ok<AGROW>
    end

    % Calculate final mean and SEM
    finalMeanPSD = mean(aggregatedMeanPSD, 1);
    finalSEMPSD = sqrt(sum(aggregatedSEMPSD.^2, 1)) / sqrt(size(aggregatedSEMPSD, 1));

    % Generate final plot
    figure;
    errorbar(combinedTimeBins, finalMeanPSD, finalSEMPSD, '-o', 'LineWidth', 2);
    xlabel('Time (hours)');
    ylabel('Normalized Percent Band Power');
    title('Normalized Percent Band Power Across Animals');
    grid on;

    % Save plot as PNG
    outputPlot = strrep(outputFile, '.mat', '.png');
    saveas(gcf, outputPlot);
    fprintf('Saved final plot to %s\n', outputPlot);

    % Save combined data
    save(outputFile, 'finalMeanPSD', 'finalSEMPSD', 'combinedTimeBins', '-v7.3');
    fprintf('Saved combined data to %s\n', outputFile);

    % Close figure to free memory
    close;
end
