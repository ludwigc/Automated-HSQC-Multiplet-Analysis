function [outliers] = find_outliers (vector)
% Compute the median absolute difference
meanValue = mean(vector);
% Compute the absolute differences.  It will be a vector.
absoluteDeviation = abs(vector - meanValue);
% Compute the median of the absolute differences
mad = median(absoluteDeviation);
% Find outliers.  They're outliers if the absolute difference
% is more than some factor times the mad value.
sensitivityFactor = 3.0; %1.5;
thresholdValue = sensitivityFactor * mad;
outlierIndexes = abs(absoluteDeviation) > thresholdValue;
% Extract outlier values:
outliers = vector(outlierIndexes);
end