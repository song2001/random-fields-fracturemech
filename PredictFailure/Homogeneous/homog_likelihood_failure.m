% Vincente Pericoli
% UC Davis
% 7 Jan 2016
%
% this function computes the likelihood of observing the failure data set 
% (samples struct), under the assumptions of deterministic l* with 
% homogeneous material properties defined by distType and distParams.
%

function [lkhood, samples, lkhoods] = homog_likelihood_failure ...
                    (samples, lstarIndex, distType, distParams, saveStruct)
%HOMOG_LIKELIHOOD_FAILURE
% INPUTS--
%   samples    = struct of the samples (the deterministic VGPy_database)
%   lstarIndex = int l* index, where lstar(lstarIndex) = l*
%   distType   = string name of the failure PDF type (e.g. 'Normal',
%                'Lognormal', etc.)
%   distParams = vector of the two PDF parameters [mean, stdev]
%   saveStruct = optional boolean to toggle saving info to the samples
%                struct (default = True)
% OUTPUTS--
%   lkhood  = float product of all lkhoods (i.e. the likelihood of
%             observing failure of the entire set)
%   lkhoods = struct of the sample names and their likelihoods
%

% set saveStruct
if nargin < 5, saveStruct = true; end

% add path for lagrange_diff()
addpath('..');

% want unique name for these results (for this l* index)
str_res = sprintf('results%i',lstarIndex);

% obtain the names & number of samples
sampleNames = fieldnames(samples);
numSam = length(sampleNames);

% obtain material type
material = samples.(sampleNames{1}).material;

% preallocate lkhoods vector
lkhoods = [];

% loop through every sample
for s = 1:numSam
    % check the material type
    if ~strcmpi(samples.(sampleNames{s}).material, material)
        % halt execution if they vary.
        error('Multiple Material Types Defined!');
    end
    
    % set this sample's load history for readability
    loadHist = samples.(sampleNames{s}).loadHist;
    
    % set this sample's VGI for readability
    VGI = samples.(sampleNames{s}).VGI;
    
    % obtain the size of the VGI array.
    % nrow corresponds to the number of history points.
    % ncol corresponds to the number of locations (deterministic l*).
    % only samples with strong gradients have more than one location.
    [~, ncol] = size(VGI);
    
    if ncol ~= 1
        % the l* location is treated as deterministic, so we want to
        % calculate the likelihood based on a pre-set deterministic l*, so
        % only include the VGIs corresponding to that l* location.
        VGI = VGI(:,lstarIndex);
    end
    
    % for the deterministic critical location, compute the failure CDF 
    % using the defined distType and distParams. Although this is the 
    % failure CDF for VGI, it is also the failure CDF for the observed 
    % loading, because they are related one-to-one.
    failCDF = cdf(distType, VGI, distParams(1), distParams(2));
    
    % now we have a (numerical) failure CDF for the observed loading, based 
    % on (1) the assumed VGI failure distType and distParams, and (2) the 
    % assumption that failure occurs at a specific critical location 
    % (deterministic). but what we want is the likelihood, so compute the 
    % numerical derivative of this CDF at the point of observed ("actual") 
    % failure.
    
    % assign the observed failure indices for this sample
    fi = samples.(sampleNames{s}).failureIndex;
    
    % numerically differentiate CDF using a 2nd order Lagrange 
    % interpolating polynomial, since the data is nonequispaced
    failPDF = mderiv_fornberg(1, loadHist, failCDF);
    
    % save the lkhoods to the sample struct, if requested
    if saveStruct
        samples.(sampleNames{s}).(str_res).lkhoods    = failPDF(fi);
        samples.(sampleNames{s}).(str_res).failCDF    = failCDF;
        samples.(sampleNames{s}).(str_res).failPDF    = failPDF;
        samples.(sampleNames{s}).(str_res).distType   = distType;
        samples.(sampleNames{s}).(str_res).distParams = distParams;
    end
    % append sample lkhood to lkhoods
    lkhoods((end+1):(end+length(failPDF(fi)))) = failPDF(fi)';

end

% likelihood of set is product of all likelihoods
lkhood = prod(lkhoods);

return;
end
