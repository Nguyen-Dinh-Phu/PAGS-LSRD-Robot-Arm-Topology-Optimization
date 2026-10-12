
%---------------------- ElementToNodalSensitivity.m -----------------------%
% Project element sensitivities onto Q4 nodes with optional spatial filtering.
% Input : TDN_elem   - Element sensitivity field
%         nelx, nely - Number of elements along x and y
%         Rmin       - Filter radius in element-grid units (default: 0)
% Output: TDN_nodal  - Projected nodal sensitivity vector
%-------------------------------------------------------------------------%

function TDN_nodal = ElementToNodalSensitivity(TDN_elem,nelx,nely,Rmin)
if nargin < 4, Rmin = 0; end           % Disable filtering by default
TDN_elem = reshape(TDN_elem,nely,nelx);% Element sensitivity matrix
H0 = ones(2);
S = conv2(TDN_elem,H0,'full');         % Sum adjacent element sensitivities
N = conv2(ones(nely,nelx),H0,'full');  % Number of adjacent elements
S = S./N;                              % Arithmetic nodal projection
if Rmin > 0
    r = ceil(Rmin);
    [dx,dy] = meshgrid(-r:r,-r:r);
    H = max(0,Rmin-hypot(dx,dy));      % Cone-shaped filter weights
    S = conv2(S,H,'same')./ ...
        conv2(ones(size(S)),H,'same'); % Boundary-normalized filtering
end
TDN_nodal = S(:);                      % Nodal sensitivity vector
%-------------------------------------------------------------------------%
