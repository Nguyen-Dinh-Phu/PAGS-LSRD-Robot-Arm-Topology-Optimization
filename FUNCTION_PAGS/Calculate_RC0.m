
%---------------------------- Calculate_RC0.m ----------------------------%
% Compute the reference reaction center from support conditions.
% Continuous supports use the geometric centroid; otherwise, the
% reaction-magnitude-weighted centroid is evaluated.
% Input : fem          - Finite element model and boundary conditions
%         K            - Assembled global stiffness matrix
%         U            - Global displacement vector
% Output: x_RC0, y_RC0 - Reference reaction center coordinates
%-------------------------------------------------------------------------%

function [x_RC0,y_RC0] = Calculate_RC0(fem,K,U)
fixedDofs = fem.bc.fixedDofs;
fixedNodes = unique(floor((fixedDofs+1)/2));  % Constrained node indices
coord = fem.mesh.node;
fixedCoord = coord(fixedNodes,:);             % Constrained node coordinates
nelx = fem.domain.nelx;
nely = fem.domain.nely;

isClamped = numel(fixedNodes) >= min(nelx,nely)/2 && ...
    (var(fixedCoord(:,1)) < 1e-5 || ...
     var(fixedCoord(:,2)) < 1e-5);            % Continuous support criterion
if isClamped
    RC = mean(fixedCoord,1);                  % Geometric support centroid
else
    R = K*U - fem.bc.forceVector;             % Global nodal reaction forces
    R(fem.bc.freeDofs) = 0;                   % Retain constrained reactions
    Rmag = hypot(R(1:2:end),R(2:2:end));      % Reaction magnitude per node
    Rtot = sum(Rmag);
    if Rtot > 1e-6
        RC = (Rmag'*coord)/Rtot;              % Reaction-weighted centroid
    else
        RC = [fem.domain.lx,fem.domain.ly]/2; % Fallback reference center
    end
end
x_RC0 = RC(1); y_RC0 = RC(2);
%-------------------------------------------------------------------------%
