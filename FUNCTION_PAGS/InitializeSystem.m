
%--------------------------- InitializeSystem.m --------------------------%
% Initialize orthotropic Q4 finite element and reaction-diffusion systems.
% Input : fem - Finite element model and mesh geometry
%         opt - Material properties and optimization parameters
% Output: SYS - FEM matrices and level-set operators
%         phi - Initial nodal level-set field
%         str - Initial element density field
%-------------------------------------------------------------------------%

function [SYS,phi,str] = InitializeSystem(fem,opt)
nelx = fem.domain.nelx; nely = fem.domain.nely;
lx = fem.domain.lx; ly = fem.domain.ly;
dx = lx/nelx; dy = ly/nely;                   % Element dimensions
nelem = nelx*nely;                            % Number of Q4 elements
nnodes = (nelx+1)*(nely+1);                   % Number of mesh nodes

%----------------------------- ORTHOTROPIC MATERIAL
E11   = opt.E11; E22 = opt.E22;
G12   = opt.G12; nu12 = opt.nu12;
nu21  = nu12*E22/E11;                       % Reciprocal Poisson's ratio
denom = 1-nu12*nu21;
if denom <= 0
    error('InitializeSystem:InvalidMaterial', ...
          'Invalid orthotropic material properties.');
end

Q11 = E11/denom;
Q22 = E22/denom;
Q12 = nu12*E22/denom;
Q66 = G12;
SYS.D0 = [ Q11 Q12 0;
           Q12 Q22 0; 
           0   0   Q66];                    % Local plane-stress stiffness

%----------------------------- MATERIAL ORIENTATION
theta  = opt.theta;
m      = cos(theta); n = sin(theta);
m2     = m^2; n2 = n^2;
m3     = m^3; n3 = n^3;
m4     = m^4; n4 = n^4;
Qbar11 = Q11*m4+2*(Q12+2*Q66)*m2*n2+Q22*n4;
Qbar22 = Q11*n4+2*(Q12+2*Q66)*m2*n2+Q22*m4;
Qbar12 = (Q11+Q22-4*Q66)*m2*n2+Q12*(m4+n4);
Qbar16 = (Q11-Q12-2*Q66)*m3*n-(Q22-Q12-2*Q66)*m*n3;
Qbar26 = (Q11-Q12-2*Q66)*m*n3-(Q22-Q12-2*Q66)*m3*n;
Qbar66 = (Q11+Q22-2*Q12-2*Q66)*m2*n2+Q66*(m4+n4);

SYS.D = [Qbar11 Qbar12 Qbar16; 
         Qbar12 Qbar22 Qbar26; 
         Qbar16 Qbar26 Qbar66];              % Transformed stiffness matrix
SYS.D = 0.5*(SYS.D+SYS.D');

%----------------------------- ELEMENT MATRICES
a     = dx/2; b = dy/2;
coord = [-a -b; a -b; a b; -a b];            % Q4 reference coordinates
g     = 1/sqrt(3);
gaussPts = [-g -g; g -g; g g; -g g];         % 2x2 Gauss quadrature
weights  = ones(4,1);
SYS.KE   = zeros(8,8);                       % Elastic element stiffness
Ke_diff  = zeros(4,4);                       % Diffusion element stiffness
for q = 1:4
    xi = gaussPts(q,1); eta = gaussPts(q,2);
    dN_dxi  = 0.25*[-(1-eta), (1-eta),(1+eta),-(1+eta)];
    dN_deta = 0.25*[-(1-xi),-(1+xi),(1+xi), (1-xi)];
    dN_nat  = [dN_dxi;dN_deta];               % Natural shape gradients
    J    = dN_nat*coord;                      % Isoparametric Jacobian
    detJ = det(J);

    if detJ <= 0
        error('InitializeSystem:InvalidJacobian', ...
              'Non-positive element Jacobian.');
    end

    gradN = J\dN_nat;                        % Physical shape gradients
    B     = zeros(3,8);
    B(1,1:2:end) = gradN(1,:);
    B(2,2:2:end) = gradN(2,:);
    B(3,1:2:end) = gradN(2,:);
    B(3,2:2:end) = gradN(1,:);               % Plane-stress strain matrix
    SYS.KE = SYS.KE+B'*SYS.D*B*detJ*weights(q);
    Ke_diff = Ke_diff+gradN'*gradN*detJ*weights(q);
end
SYS.KE = 0.5*(SYS.KE+SYS.KE');               % Symmetric element stiffness

Me = (dx*dy/36)*[4 2 1 2; 
                 2 4 2 1; 
                 1 2 4 2; 
                 2 1 2 4];                   % Consistent Q4 mass matrix

%----------------------------- FEM CONNECTIVITY
nodenrs = reshape(1:nnodes,nely+1,nelx+1);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelem,1);
offset  = [0 1 2*nely+2 2*nely+3 2*nely 2*nely+1 -2 -1];
SYS.edofMat = repmat(edofVec,1,8)+ ...
              repmat(offset,nelem,1);        % Displacement DOF connectivity
SYS.iK = reshape(kron(SYS.edofMat,ones(8,1))',64*nelem,1);
SYS.jK = reshape(kron(SYS.edofMat,ones(1,8))',64*nelem,1);

%----------------------------- LEVEL-SET CONNECTIVITY
edofVec_phi     = reshape(nodenrs(1:end-1,1:end-1),nelem,1);
offset_phi      = [0 nely+1 nely+2 1];
SYS.edofMat_phi = repmat(edofVec_phi,1,4)+ ...
                  repmat(offset_phi,nelem,1); % Scalar Q4 connectivity
iN = reshape(kron(SYS.edofMat_phi,ones(4,1))',16*nelem,1);
jN = reshape(kron(SYS.edofMat_phi,ones(1,4))',16*nelem,1);

%----------------------------- REACTION-DIFFUSION OPERATORS
SYS.NN    = sparse(iN,jN,repmat(Me(:),nelem,1),nnodes,nnodes);
SYS.Kdiff = sparse(iN,jN,repmat(Ke_diff(:),nelem,1), ...
                   nnodes,nnodes);             % Global diffusion stiffness
SYS.NN    = 0.5*(SYS.NN+SYS.NN');              % Consistent mass matrix
SYS.Kdiff = 0.5*(SYS.Kdiff+SYS.Kdiff');
SYS.T     = SYS.NN/opt.dt+opt.tau*SYS.Kdiff;   % Implicit LS-RD operator
SYS.fixeddofsphi = [];
SYS.freedofsphi  = (1:nnodes)';                % Free level-set DOFs

%----------------------------- INITIAL DESIGN
phi = ones(nnodes,1);                         % Initial solid level-set field
str = ones(nely,nelx);                        % Initial element densities
if isfield(fem,'bc')
    if isfield(fem.bc,'passive_solid')
        str(fem.bc.passive_solid) = 1;        % Prescribed solid regions
    end
    if isfield(fem.bc,'passive_void')
        str(fem.bc.passive_void) = 0;         % Prescribed void regions
    end
end

%----------------------------- SYSTEM INFORMATION
SYS.nelem  = nelem;
SYS.nnodes = nnodes;
SYS.dx     = dx; SYS.dy = dy;
SYS.theta  = theta;
SYS.nu21   = nu21;
SYS.material = struct( ...
    'E11',E11, ...
    'E22',E22, ...
    'G12',G12, ...
    'nu12',nu12, ...
    'nu21',nu21);
%-------------------------------------------------------------------------%
