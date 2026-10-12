
%--------------------------- ApplyPassivePhi.m ---------------------------%
% Enforce prescribed level-set values at nodes of passive Q4 elements.
% Input : phi        - Nodal level-set field
%         elemID     - Indices of passive elements
%         value      - Prescribed level-set value
%         nelx, nely - Number of elements along x and y
% Output: phi        - Updated nodal level-set field
%-------------------------------------------------------------------------%

function phi = ApplyPassivePhi(phi,elemID,value,nelx,nely)
[iy,ix] = ind2sub([nely,nelx],elemID(:));    % Element subscripts
n = (ix-1)*(nely+1) + iy;                    % Lower-left node indices
nodes = [n;n+1;n+nely+1;n+nely+2];           % Four corner nodes
phi(unique(nodes)) = value;                  % Enforce passive constraints
%-------------------------------------------------------------------------%
