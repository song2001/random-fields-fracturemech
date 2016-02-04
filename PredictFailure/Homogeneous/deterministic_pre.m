% Vincente Pericoli
% UC Davis
% 21 Jan 2016
%
% preprocess the "deterministic/homogeneous" results. returns a struct 
% containing only the requested material
%

function [mySamples, lstars] = deterministic_pre(samples, material)
%PREPROCESS - pre-process the deterministic results
%
%Inputs-
%   samples  = struct of the samples (the deterministic VGPy_database)
%   material = string name of the requested material
%

% set user material definition into cell for convenience
if strcmpi(material,'AP70HP') || strcmpi(material,'HPS70W')
    material = {'AP70HP','HPS70W'};
else
    material = {material};
end

% obtain the names & number of samples
sNames = fieldnames(samples);
nums = length(sNames);

% loop through all samples
for s = 1:nums
    
    % if the sample is BB or BH, a single VGI (max) needs to be chosen for
    % each time history (frame)
    if strncmpi(sNames{s},'BB_',3) || strncmpi(sNames{s},'BH_',3)
        % this sample is a BB or BH--
        VGI = samples.(sNames{s}).VGI;
        samples.(sNames{s}).VGI = max(VGI,[],2);
    end
    
    % also return sample struct containing ONLY the requested material type
    if any( strcmpi(samples.(sNames{s}).material, material) )
        % this is the requested material, save to mySamples struct
        mySamples.(sNames{s}) = samples.(sNames{s});
    end
    
    % we also want find out what the possible l*'s of this material are
    if strncmpi(sNames{s},'CT_',3)
        lstars = samples.(sNames{s}).lstars;
    end
end


return;
end