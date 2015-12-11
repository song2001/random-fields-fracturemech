% Vincente Pericoli
% UC Davis
% random-fields-fracturemech
% 8 Dec 2015
%
% Contains function gauss() and subfunction gauss_rule()
% used to return the gauss weights and coordinates for quadrature
%

function [gw, gpt] = gauss(ndim, nippe)
%GAUSS
% Calculates the gauss weights and coordinates for all integration points
% for the parent domain (i.e. contained in [-1,1]^ndim), assuming the 
% ABAQUS integration point numbering convention (see the ABAQUS 
% documentation for more info)
%
%Inputs -
%   ndim  : integer number of dimensions (2 or 3)
%   nippe : integer number of integration points per element
%
%Outputs -
%   gw  : vector of gauss weights, gw(i) is the gauss weight associated 
%         with integration point i
%   gpt : array of the gauss points (i.e. parent domain coordinates), 
%         gpt(i,:) are the coordinates associated with integration point i
%

% calculate the number of sampling points, AKA Gauss "rule"
nsp = nippe^(1/ndim);

% retrieve the the gauss points and weights associated with this rule
[p, w] = gauss_rule(nsp);

% preallocate/initialize
gw  = zeros(nippe,1);    % gauss weights for the int.pts.
gpt = zeros(nippe,ndim); % (gauss) point locations of the int.pts.
n   = 1;                 % dummy index variable

if ndim == 2
    % if 2-dimensional element, we only have two relevant coordinates, so
    % only two for loops
    for j = 1:nsp
        % corresponds to y-coordinate
        for i = 1:nsp
            % corresponds to x-coordinate
            gw(n)    = w(i)*w(j);
            gpt(n,:) = [p(i), p(j)];
            n = n + 1;
        end
    end
    
elseif ndim == 3
    % if 3-dimensional element, need three for loops
    for k = 1:nsp
        % corresponds to z-coord
        for j = 1:nsp
            % corresponds to y-coord
            for i = 1:nsp
                % corresponds to x-coord
                gw(n)    = w(i)*w(j)*w(k);
                gpt(n,:) = [p(i), p(j), p(k)];
                n = n + 1;
            end
        end
    end
else
    error('Illegal dimension');
end
    
end


function [p, w] = gauss_rule(nsp)
%GAUSS_RULE: Gauss-Legendre quadrature rule over the interval [-1,1]
% Given the number of sampling points (nsp <= 3), returns a vector of the
% Gauss points and a corresponding vector of the Gauss weights

% preallocate
p = zeros(nsp,1);
w = zeros(nsp,1);

% retrieve points and weights
switch nsp
    case 1
        % points
        p(1) =   0.0;
        % weights
        w(1) =   2.0;
    case 2
        % points
        p(1) = - 0.577350269189625764509148780502;
        p(2) =   0.577350269189625764509148780502;
        % weights
        w(1) =   1.0; 
        w(2) =   1.0;
    case 3
        % points
        p(1) = - 0.774596669241483377035853079956; 
        p(2) =   0.0;
        p(3) =   0.774596669241483377035853079956;
        % weights
        w(1) =   0.5555555555555555555555555555565;
        w(2) =   0.8888888888888888888888888888889;
        w(3) =   0.5555555555555555555555555555565;
    otherwise
        error('Illegal number of sampling points!')
end

end