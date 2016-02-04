% Vincente Pericoli
% UC Davis
%
% this set of functions uses optimization techniques to choose the
% optimal distparams for achieving maximum likelihood, given a specific
% distType and deterministic l*.
%

function [distParams] = homog_optim_likelihood ...
                                   (samples, lstarIndex, distType, params0)
%DETERMINISTIC_OPTIM_LIKELIHOOD
% Uses the Nelder-Mead Simplex Method to find distParams which achieve the
% highest possible likelihood. This inherrently assumes that the space is
% convex or well-behaved.
%
% Inputs--
%   samples
%   lstarIndex
%   distType
%   params0 = vector initial guess of distParams (default = [0,1])
% Outputs--
%   distParams
%

% use persistent optimization options (much faster if this is in a loop)
persistent OPTIONS
if isempty(OPTIONS)
    OPTIONS = optimset(@fminsearch);
    OPTIONS.FunValCheck = 'on';
    OPTIONS.MaxFunEvals = 1000;
    OPTIONS.MaxIter     = 1000;
    OPTIONS.PlotFcns    = @optimplotfval;
end

% assign params0 if not already defined (assumes 2-parameter distType)
if nargin < 4, params0 = [0,1]; end

% perform optimization using Nelder-Mead Simplex Method
opt = fminsearch(@(X) objective_fun(X,samples,lstarIndex,distType), ...
                            params0, OPTIONS);

distParams = opt;

return;
end

function [lkhood] = objective_fun(distParams,samples,lstarIndex,distType)
% objective function

% given the distParams, calculate likelihood
[lkhood, ~, ~] = homog_likelihood_failure ... 
                               (samples, lstarIndex, distType, distParams);
                           
% we want to find the maximum likelihood, but optimization is an argmin().
% so, perform optimization on the negative. this is equivalent to finding
% the argmax(), assuming that the space is convex.
lkhood = -lkhood;

return;
end