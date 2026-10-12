
%--------------------------- PhiToDensity.m ------------------------------%
% Map nodal level-set values to element densities using a regularized
% Heaviside function.
% Input : phi        - Nodal level-set field
%         nelx, nely - Number of elements along x and y
%         delta      - Regularization half-width
% Output: str        - Element density field in [0,1]
%-------------------------------------------------------------------------%

function str = PhiToDensity(phi,nelx,nely,delta)
phin = reshape(phi,nely+1,nelx+1);       % Nodal level-set field
phie = 0.25*(phin(1:end-1,1:end-1) + ...
             phin(2:end,1:end-1) + ...
             phin(1:end-1,2:end) + ...
             phin(2:end,2:end));         % Element-averaged level-set values
str = zeros(size(phie));                 % Initialize void density
str(phie > delta) = 1;                   % Solid region

idx = abs(phie) <= delta;                % Regularized transition region
xi = phie(idx)/delta;
str(idx) = 0.5*(1 + xi + sin(pi*xi)/pi); % Smooth Heaviside interpolation
%-------------------------------------------------------------------------%
