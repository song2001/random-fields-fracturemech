% Vincente Pericoli
% UC Davis
%
% driver to plot the likelihood parameter space given a set of means and
% standard deviations
%

% clear workspace
clear; clc; close all;

% change dirs, load data, change back to pwd
fdir = pwd;
addpath('..');
cd(myPaths('VGPy_Databases'));
samples = load('Deterministic.mat');
cd(fdir); clear fdir;

% perform any necessary preprocessing
[~, samplesAP50, samplesAP70HP] = deterministic_pre(samples);

% choose the material
samples = samplesAP50;

% define the desired distribution
distType = 'Lognormal';

% brute force means, stdevs to plot over
means  = 0.1:0.1:10;
stdevs = 0.1:0.1:5;

for lstarIndex = 1:3
    % use brute-force to determine appropriate distParams
	[map, ~, mean, stdev] = deterministic_brute_likelihood ...
                        (samplesAP50, lstarIndex, distType, means, stdevs);
    % map the parameter space
	figure;
	surf(stdevs,means,map);
    title(sprintf('l* index = %i',lstarIndex));
    xlabel('stdev'); ylabel('mean'); zlabel('likelihood');
end