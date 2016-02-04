% Vincente Pericoli
% UC Davis
%
% given vectors of means and stdevs, provides a brute-force map of the
% likelihood parameter space. Also will provide the best likelihood
% achieved from the provided means and stdevs. Very slow.
%

function [map, bestval, bestmean, beststdev] = ...
            homog_param_space(samples, lstarIndex, distType, means, stdevs)
%HOMOG_PARAM_SPACE
%

bestmean  = -1;
beststdev = -1;
bestval = -inf;

map = zeros(length(means),length(stdevs));

for m = 1:length(means)
    for s = 1:length(stdevs)
        % set params
        distParams(1) = means(m);
        distParams(2) = stdevs(s);
        % get lkhood associated with params
        [lkhood, ~,~] = deterministic_likelihood_failure ...
                               (samples, lstarIndex, distType, distParams);
        % save to map
        map(m,s) = lkhood;       
        % if this lkhood better than current bestval, replace.
        if lkhood > bestval
            bestval   = lkhood;
            bestmean  = means(m);
            beststdev = stdevs(s);
        end
    end
end

return;
end