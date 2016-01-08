% Vincente Pericoli
% UC Davis
% 7 Jan 2016


function [lkhood, lkhoods] = deterministic_likelihood_failure ... 
                                (samples, LstarIndex, distType, distParams)
%DETERMINISTIC_LIKELIHOOD_FAILURE
% INPUTS--
%   samples    = struct of the samples (the deterministic VGPy_database)
%   LstarIndex = int l* index, where Lstar(LstarIndex) = l*
%   distType   = string name of the failure PDF type (e.g. 'Normal',
%                'Lognormal', etc.)
%   distParams = vector of the two PDF parameters [mean, stdev]
%
% OUTPUTS--
%   lkhood  = float product of all lkhoods (i.e. the likelihood of
%             observing failure of the entire set)
%   lkhoods = vector of the likelihoods for each sample
%

% obtain the names of the samples
sampleNames = fieldnames(samples);
numSam = length(sampleNames);

% preallocate lkhoods
lkhoods = NaN(1,numSam);

% loop through every sample
for s = 1:numSam
    
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
        VGI = VGI(:,LstarIndex);
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
    
    % assign the failure index for this sample
    fi = samples.(sampleNames{s}).failureIndex;
    
    % numerically differentiate using a 2nd order Lagrange interpolating
    % polynomial, since the data is nonequispaced
    w1 = loadHist(fi) - loadHist(fi+1);
    w1 = w1/( (loadHist(fi-1) - loadHist(fi)) * ...
              (loadHist(fi-1) - loadHist(fi+1)) );
	w2 = 2*loadHist(f1) - loadHist(fi-1) - loadHist(fi+1);
    w2 = w2/( (loadHist(fi) - loadHist(fi-1)) * ...
              (loadHist(fi) - loadHist(fi+1)) );
	w3 = loadHist(fi) - loadHist(fi-1);
    w3 = w3/( (loadHist(fi+1) - loadHist(fi-1)) * ...
              (loadHist(fi+1) - loadHist(fi)) );
	
    lkhoods(s) = failCDF(fi-1)*w1 + failCDF(fi)*w2 + failCDF(fi+1)*w3;
    
end

lkhood = prod(lkhoods);
end