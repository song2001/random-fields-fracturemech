% Vincente Pericoli
% UC Davis
%
% homogeneous material driver
%

% clear workspace
clear; clc; close all;

% change dirs, load data, change back to pwd
fdir = pwd;
addpath('..');
cd(myPaths('VGPy_Databases'));
samples = load('Deterministic.mat');
cd(fdir);

% user inputs
distType = 'Lognormal';
saveflag = false;
material = 'AP50';

% perform any necessary preprocessing
[samples, lstars] = deterministic_pre(samples, material);

% initialize
bestlkhood = -inf;
bestlstar  = -1;

for lstarIndex = 1:length(lstars)
    % use optimization to determine appropriate distParams
    distParams = homog_optim_likelihood(samples, lstarIndex, distType);
    
    % determine the likelihoods for this set, using distParams
    % this also returns a samples with a unique results field corresponding
    % to this particular LstarIndex
    [lkhood, samples, ~] = homog_likelihood_failure ...
                               (samples, lstarIndex, distType, distParams);
    % keep track of the "best" candidate l*
    if lkhood > bestlkhood
        bestlstar  = lstars(lstarIndex);
        bestlkhood = lkhood;
    end
    
    % plot the CDFs for each sample
    %plot_cdfs(samples,LstarIndex,lstars,saveflag);
    %plot_pdfs(samples,LstarIndex,lstars,saveflag);
    close all;
    
    % save these results to a struct to output
    lsn = sprintf('lstar%i',lstarIndex);
    results.(material).(lsn).distParams = distParams;
    results.(material).(lsn).lkhood     = lkhood;
    results.(material).(lsn).lstar      = lstars(lstarIndex);
end

% save final samples & best info to results struct
results.(material).samples     = samples;
results.(material).best.lkhood = bestlkhood;
results.(material).best.lstar  = bestlstar;

%cd(myPaths('save-results-homogeneous'));
%save('results_optim','results');
%cd(fdir);