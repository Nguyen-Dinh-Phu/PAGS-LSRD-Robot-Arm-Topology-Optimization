%------------------------------- ApplyBC.m -------------------------------%
% Define boundary conditions and passive regions for 2D FEM models.
% Cases : '1' MBB beam, '2' Cantilever beam, '3' Robotic arm
% Input : fem    - Finite element model
%         typeBC - Boundary condition identifier
%         F      - Applied force magnitude
% Output: fem    - Updated supports, loads, and passive regions
%-------------------------------------------------------------------------%

function fem = ApplyBC(fem,typeBC,F)
nelx   = fem.domain.nelx; nely = fem.domain.nely;
Lx     = fem.domain.lx; Ly = fem.domain.ly;
dx     = Lx/nelx; dy = Ly/nely;
nNodes = (nelx+1)*(nely+1);
nDofs  = 2*nNodes;
forceVector = sparse(nDofs,1);              % Global nodal load vector
fixedDofs   = [];
fem.bc.passive_void  = [];                  % Prescribed void elements
fem.bc.passive_solid = [];                  % Prescribed solid elements

switch lower(char(typeBC))
    case '1' % MBB BEAM
        leftNodes = 1:nely+1;
        rightNode = nelx*(nely+1)+1;
        fixedDofs = [2*leftNodes-1,2*rightNode]; % Symmetry and roller supports
        forceVector(2) = -F;                     % Downward point load

    case '2' % CANTILEVER BEAM
        leftNodes  = 1:nely+1;
        rightNodes = nelx*(nely+1)+1:(nelx+1)*(nely+1);
        loadNode   = rightNodes(round(numel(rightNodes)/2));
        fixedDofs  = [2*leftNodes-1,2*leftNodes]; % Clamped left edge
        forceVector(2*loadNode) = -F;             % Downward point load

    case '3' % ROBOTIC ARM
        r_ring_A = 15; r_hole_A = 8;
        r_ring_B = 13; r_hole_B = 6;
        xcA = r_ring_A;       ycA = Ly/2;        % Left joint center
        xcB = Lx-r_ring_B-2;  ycB = Ly/2;        % Right joint center
        [Xe,Ye] = meshgrid(dx/2:dx:Lx-dx/2, ...
                           dy/2:dy:Ly-dy/2);     % Element centroids
        dAe = hypot(Xe-xcA,Ye-ycA);
        dBe = hypot(Xe-xcB,Ye-ycB);              % Distances to joint centers
        outerVoid = (Xe<xcA & dAe>r_ring_A);
        innerVoid = (dAe<r_hole_A | dBe<r_hole_B);
        fem.bc.passive_void  = find(outerVoid | innerVoid);
        fem.bc.passive_solid = find( (dAe>=r_hole_A & dAe<=r_ring_A) | ...
                                         (dBe>=r_hole_B & dBe<=r_ring_B));
        [Xn,Yn] = meshgrid(linspace(0,Lx,nelx+1), ...
                           linspace(0,Ly,nely+1));  % Nodal coordinates
        nodeID = reshape(1:nNodes,nely+1,nelx+1);
        tol = max(dx,dy);                           % Rim selection tolerance
        dAn = hypot(Xn-xcA,Yn-ycA);
        supportNodes = nodeID(abs(dAn-r_hole_A)<=tol);
        fixedDofs    = [2*supportNodes(:)-1; ...
                        2*supportNodes(:)];             % Fixed inner rim of joint A
        dBn = hypot(Xn-xcB,Yn-ycB);
        rim = abs(dBn-r_hole_B)<=tol;
        loadNodes = nodeID(rim & Yn<=ycB);              % Lower inner rim of joint B
        forceVector(2*loadNodes) = -F/numel(loadNodes); % Distributed vertical load
end
fem.bc.fixedDofs = unique(fixedDofs);           
fem.bc.freeDofs = setdiff(1:nDofs,fem.bc.fixedDofs);
fem.bc.forceVector = forceVector;               
%-------------------------------------------------------------------------%
