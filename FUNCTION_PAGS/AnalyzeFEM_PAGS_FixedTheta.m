
%---------------------- AnalyzeFEM_PAGS_FixedTheta.m ----------------------%
% Evaluate compliance, CG-RC alignment, and design sensitivities.
% Input : str      - Element density field
%         fem, SYS - Finite element model and system matrices
%         OPT      - Material and optimization parameters
%         MeshInfo - Element areas and centroid coordinates
%         C0, d02  - Reference compliance and squared CG-RC distance
% Output: objective, TDN_total - Objective and descent sensitivity
%         C_val, d2_val        - Compliance and squared CG-RC distance
%         x_CG, y_CG, Iz_val   - Geometric response quantities
%-------------------------------------------------------------------------%

function [objective,TDN_total,C_val,d2_val,x_CG,y_CG,Iz_val] = ...
    AnalyzeFEM_PAGS_FixedTheta(~,str,fem,SYS,OPT,MeshInfo,C0,d02)
nelx = fem.domain.nelx;
nely = fem.domain.nely;
nelem = nelx*nely;
rho = max(0,min(1,str(:)));                   
Emin = OPT.Emin;
p = OPT.penal;
w = OPT.w;

%----------------------------- MATERIAL INTERPOLATION
E = Emin + (1-Emin)*rho.^p;                   % Ersatz stiffness factors

%----------------------------- FINITE ELEMENT ANALYSIS
sK = reshape(SYS.KE(:)*E',64*nelem,1);        % Element stiffness contributions
K = sparse(SYS.iK,SYS.jK,sK);                 % Global stiffness matrix
K = 0.5*(K+K');                               % Enforce matrix symmetry
freeDofs = fem.bc.freeDofs;
U = zeros(2*(nelx+1)*(nely+1),1);
U(freeDofs) = K(freeDofs,freeDofs) \ ...
              fem.bc.forceVector(freeDofs);   % Equilibrium displacement field
Ue = U(SYS.edofMat);                          % Element displacement vectors
ce = sum((Ue*SYS.KE).*Ue,2);                  % Element strain-energy terms
C_val = sum(E.*ce);                           % Structural compliance
rhoSafe = max(rho,1e-8);
dCdrho = -p*(1-Emin)*rhoSafe.^(p-1).*ce;      % Compliance sensitivity

%----------------------------- CENTER OF GRAVITY
Area = MeshInfo.Area(:);
xc = MeshInfo.Centroid(:,1);
yc = MeshInfo.Centroid(:,2);
TotalMass = sum(rho.*Area);                   % Density-weighted area measure

if TotalMass > 1e-12
    x_CG = sum(rho.*Area.*xc)/TotalMass;      % CG x-coordinate
    y_CG = sum(rho.*Area.*yc)/TotalMass;      % CG y-coordinate
else
    x_CG = fem.domain.lx/2;
    y_CG = fem.domain.ly/2;
end

%----------------------------- CG-RC ALIGNMENT
x_RC = OPT.CR(1);
y_RC = OPT.CR(2);
dx = x_CG-x_RC;
dy = y_CG-y_RC;
d2_val = dx^2+dy^2;                           % Squared CG-RC distance
if TotalMass > 1e-12
    dxCGdrho = Area.*(xc-x_CG)/TotalMass;     % CG sensitivity along x
    dyCGdrho = Area.*(yc-y_CG)/TotalMass;     % CG sensitivity along y
else
    dxCGdrho = zeros(nelem,1);
    dyCGdrho = zeros(nelem,1);
end
dd2drho = 2*dx.*dxCGdrho+2*dy.*dyCGdrho;      % Alignment sensitivity

%----------------------------- MULTI-OBJECTIVE FUNCTION
Cref = max(abs(C0),1e-12);
Dref = max(abs(d02),1e-12);
objective = (1-w)*C_val/Cref+w*d2_val/Dref;   % Normalized objective
dJdrho = (1-w)*dCdrho/Cref+ ...
                           w*dd2drho/Dref;    % Combined density sensitivity

%----------------------------- GEOMETRIC RESPONSE
r2 = (xc-x_RC).^2+(yc-y_RC).^2;
Iz_val = sum(rho.*Area.*r2);                   % Polar second-moment proxy

%----------------------------- LEVEL-SET DRIVING FIELD
TDN_total = reshape(-dJdrho,nely,nelx);        % Descent sensitivity field
%-------------------------------------------------------------------------%
