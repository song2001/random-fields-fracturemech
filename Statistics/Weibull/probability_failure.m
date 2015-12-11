% Vincente Pericoli
% UC Davis
% random-fields-fracturemech
% 8 Dec 2015
%
% Weibull failure CDF
%

function pfail = probability_failure ...
                    (VGI_ELEM_IP, lmtype, elemConnect, nodesCoords, params)
%PROBABILITY_FAILURE
% This function calculates the Weibull failure CDF for a given specimen.
%
%Inputs -
%
%   VGI_ELEM_IP : rank-3 array of the VGI at integration points (IPs) for 
%                 each element. access should correspond to (frame,IP,elem)
%
%   lmtype      : string of the ABAQUS element type
%
%   elemConnect : rank-2 array of the element connectivity. access should
%                 correspond to (elem,node), where the nodes are numbered
%                 with ABAQUS convention
%
%   nodesCoords : rank-2 array of the nodal coordinates. access should
%                 correspond to (node,dimension)
%
%   params      : vector of Weibull fitting parameters. Can be of length 3
%                 or 4. Valid inputs can be one of the following:
%                 [Weibull modulus, VGI threshold, (V0)*(VGI0)^m]
%                 [Weibull modulus, VGI threshold, VGI0, V0]
%Outputs -
%
%   pfail : vector of the failure probability at each frame
%


%
% add paths
%
addpath('..')
addpath(myPaths('fem-interp'));

%
% Determine the scope of the problem
%

% determine the element type and problem dimensionality
if strncmpi(lmtype,'CAX4',4)
    % this is a 2D axisymmetric QUAD4 element
    ndim   = 2;
    axisym = true;
    lmtype = 'QUAD4';
elseif strncmpi(lmtype,'CAX8',4)
    % this is a 2D axisymmetric QUAD8 element
    ndim   = 2;
    axisym = true;
    lmtype = 'QUAD8';
elseif strncmpi(lmtype,'CPE4',4)
    % this is a 2D plane-strain QUAD4 element
    ndim   = 2;
    axisym = false;
    lmtype = 'QUAD4';
elseif strncmpi(lmtype,'CPE8',4)
    % this is a 2D plane-strain QUAD8 element
    ndim   = 2;
    axisym = false;
    lmtype = 'QUAD8';
elseif strncmpi(lmtype,'C3D8',4)
    % this is a 3D BRK8 element
    ndim   = 3;
    axisym = false;
    lmtype = 'BRK8';
elseif strncmpi(lmtype,'C3D20',5)
    % this is a 3D BRK20 element
    ndim   = 3;
    axisym = false;
    lmtype = 'BRK20';
else
    error('unknown element type');
end

% obtain # of history frames, # of IPs per elem, # of total elems
[nHist,nippe,nele] = size(VGI_ELEM_IP);

% check nele consistency, and obtain the number of nodes per elem
[dummy,nnpe] = size(elemConnect);
if dummy ~= nele
    error('input dimensions do not agree')
end


%
% Perform some necessary pre-processing tasks
%

% calculate the relevant gauss quantities for this type of element, for the
% parent interval [-1,1]
% W(i)    = weight for integration point i
% XI(i,:) = parent space coordinates of integration point i
[W, XI] = gauss(ndim, nippe);

% obtain the basis function and derivative for ip parent space locations
N  = zeros(nnpe,nippe);
dN = zeros(nnpe,3,nippe);
for ip = 1:nippe
    [N(:,ip),dN(:,:,ip),~] = elementBasisFns( XI(ip,:), lmtype );
end

% preallocate integration results
I = zeros(nHist,nele);

% set the Weibull parameters for readability
if length(params) == 3
    % 3-parameter Weibull
    m      = params(1); % Weibull modulus
    VGIth  = params(2); % threshold VGI
    CombV0 = params(3); % a "combined" volume parameter: (V0)*(VGI0)^m
    
    % these parameters need to remain unset
    VGI0   = NaN;
    V0     = NaN;
elseif length(params) == 4
    % 4-parameter Weibull
    m      = params(1); % Weibull modulus
    VGIth  = params(2); % threshold VGI
    VGI0   = params(3); % characteristic VGI
    V0     = params(4); % characteristic sampling volume
    
    % these parameters need to remain unset
    CombV0 = NaN;
else
    error('Illegal size of params');
end
    

%
% Perform the integration
%

% loop over all of the loading history
for h = 1:nHist
    % loop over all elements, calculating the Weibull integration via Gauss
    % quadrature rule (which is consistent with that element) and the
    % isoparametric mapping.
    for e = 1:nele
        % for this element, loop over the integration points
        for ip = 1:nippe
            
            % perform a check on VGI threshold. per Weibull, if VGI for
            % this ip is below the threshold, it does not get counted.
            if VGI_ELEM_IP(h,ip,e) < VGIth
                % skip this ip (it's contribution will be zero)
                continue
            end
            
            % to use isoparametric map, we need to calculate the Jacobian.
            J = zeros(ndim,ndim);

            for k = 1:ndim
                for q = 1:ndim
                    J(k,q) = dN(:,k,ip)' * nodesCoords(elemConnect(e,:),q);
                end
            end
            
            % calculate the Weibull portion of the integrand for this ip
            weint = ( VGI_ELEM_IP(h,ip,e) - VGIth ) ^ m ;
            
            % per the quadrature rule, figure out this ip's contribution
            I_ip = W(ip) * weint * det(J);
            
            if axisym
                % if axisymmetric, the contribution needs to have a 2*pi*r
                I_ip = I_ip * ...
                       (2*pi)*(N(:,ip)' * nodesCoords(elemConnect(e,:),1));
            end
            
            % sum this into the other ip's for this element
            I(h,e) = I(h,e) + abs(I_ip);
        end
    end
end


%
% Calculate the Weibull failure CDF
%

if length(params) == 3
    % 3-parameter Weibull
    pfail = 1 - exp( - 1/CombV0 .* sum(I,2) );
    
elseif length(params) == 4
    % 4-parameter Weibull
    pfail = 1 - exp( - 1/V0/VGI0 .* sum(I,2) );

else
    % undefined
    pfail = NaN;
end

end
