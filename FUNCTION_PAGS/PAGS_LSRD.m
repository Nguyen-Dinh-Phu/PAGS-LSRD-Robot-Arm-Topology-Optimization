
%----------------------------- PAGS_LSRD.m -------------------------------%
% Multi-objective reaction-diffusion level-set topology optimization.
% Organizational structure inspired by PolyTop (Talischi et al.).
%
% Input : mesh    - Finite element model and mesh information
%         opt     - Material and optimization parameters
%         caseDef - Material orientation, objective weight, and output flags
% Output: result  - Optimized topology and numerical response
%-------------------------------------------------------------------------%

function result = PAGS_LSRD(mesh,opt,caseDef)
opt.theta = deg2rad(caseDef.thetaDeg);          % Material orientation [rad]
opt.w     = caseDef.w;                          % Multi-objective weight
opt.penal = caseDef.referencePenal;             % Reference interpolation exponent
fem       = mesh.fem;

[SYS,phi,str] = InitializeSystem(fem,opt);      % Initialize FEM and level-set fields
str           = enforcePassiveDensity(str,fem); % Enforce non-design regions
V0 = sum(str(:).*mesh.area);                    % Initial material area
[C0,D02,RC]   = ReferenceState(str,fem,SYS,opt,mesh);
opt.CR        = RC;                             % Reference reaction center
iterMax       = opt.Inter;                      % Number of optimization iterations
h = allocateHistory(iterMax,caseDef.recordHistory);
snapshotIter  = [1 100 150 300 400];            % Selected topology iterations
snapshotPhi   = cell(1,numel(snapshotIter));

if ~caseDef.recordSnapshots; snapshotPhi = {}; end
fprintf('theta=%g deg, w=%g, alpha=%g, delta=%g, p0=%g\n', ...
    caseDef.thetaDeg,opt.w,opt.alpha,opt.delta,caseDef.referencePenal);

%----------------------------- OPTIMIZATION LOOP
for k = 1:iterMax
    opt.penal = min(opt.penalMax, ...
                     1.0+(k-1)*opt.penalStep);     % Penalization continuation
    [J,dJ,C,D2] = ObjectiveFnc(phi,str,fem,SYS,opt,mesh,C0,D02);
    vol         = ConstraintFnc(str,mesh.area,V0); % Current volume fraction
    if caseDef.recordHistory
        h.iter(k) = k; h.J(k) = J; h.C(k) = C;
        h.D2(k) = D2; h.volume(k) = vol;
    end
    [phi,str] = UpdateScheme(phi,str,dJ,vol,k,fem,SYS,opt,mesh);
    if caseDef.recordSnapshots
        idx = find(snapshotIter==k,1);
        if ~isempty(idx)
            snapshotPhi{idx} = phi;             % Store selected topologies
        end
    end
    if mod(k,50)==0
        fprintf('Iter %3d/%d | C=%.8g | D2=%.8g | V=%.6f\n', ...
            k,iterMax,C,D2,vol);
    end
end

[Jfinal,~,Cfinal,D2final] = ObjectiveFnc(phi,str,fem,SYS,opt,mesh,C0,D02);
result = struct( ...
    'theta_deg',caseDef.thetaDeg, ...
    'w',caseDef.w, ...
    'alpha',opt.alpha, ...
    'delta',opt.delta, ...
    'C0',C0, ...
    'd02',D02, ...
    'xRC',RC(1), ...
    'yRC',RC(2), ...
    'Cfinal',Cfinal, ...
    'd2final',D2final, ...
    'Cnorm_final',Cfinal/max(C0,opt.tol), ...
    'Dnorm_final',sqrt(D2final/max(D02,opt.tol)), ...
    'Jfinal',Jfinal, ...
    'phi_final',phi, ...
    'finalPenal',opt.penal, ...
    'h_iter',h.iter, ...
    'h_obj',h.J, ...
    'h_comp',h.C, ...
    'h_d2',h.D2, ...
    'h_vol',h.volume, ...
    'snapshotIter',snapshotIter, ...
    'snapshotPhi',{snapshotPhi});
end

%----------------------------- OBJECTIVE FUNCTION
function [J,sensitivity,C,D2] = ObjectiveFnc(phi,str,fem,SYS,opt,mesh,C0,D02)
[J,sensitivity,C,D2] = AnalyzeFEM_PAGS_FixedTheta( ...
                                    phi,str,fem,SYS,opt,mesh.info,C0,D02);
end

%----------------------------- CONSTRAINT FUNCTION
function vol = ConstraintFnc(str,area,V0)
vol = sum(str(:).*area)/V0;                    % Normalized material volume
end

%----------------------------- LEVEL-SET UPDATE SCHEME
function [phi,str] = UpdateScheme(phi,str,sensitivity,vol,k,fem,SYS,opt,mesh)
phiOld = phi;                                  % Previous level-set field
nodal  = ElementToNodalSensitivity( ...
             sensitivity,mesh.nelx,mesh.nely); % Element-to-node projection
maxSens   = max(abs(nodal(:)));                % Maximum sensitivity magnitude
if maxSens > opt.tol
    nodal = nodal/maxSens;
else
    nodal = zeros(size(nodal));
end
[phiNew,~] = LVS(phi,nodal,vol,1.0,k,opt.nvol,opt.Vmaxin, ...
    opt.betaVol,opt.dVol,SYS.NN,opt.dt,SYS.T, ...
    SYS.freedofsphi,mesh.nelx,mesh.nely);      % Reaction-diffusion update
phi = min(1,max(-1, ...
    opt.alpha*phiNew+(1-opt.alpha)*phiOld));   % Relaxation and projection
if isfield(fem.bc,'passive_void')
    phi = ApplyPassivePhi(phi,fem.bc.passive_void, ...
        -1,mesh.nelx,mesh.nely);
end
if isfield(fem.bc,'passive_solid')
    phi = ApplyPassivePhi(phi,fem.bc.passive_solid, ...
        1,mesh.nelx,mesh.nely);
end
str = PhiToDensity(phi,mesh.nelx,mesh.nely,opt.delta);
str = enforcePassiveDensity(str,fem);          % Restore passive densities
end

%-----------------------------REFERENCE FEM STATE
function [C0,D02,RC] = ReferenceState(str,fem,SYS,opt,mesh)
rho = str(:);                                  % Initial element densities
E   = opt.Emin+(1-opt.Emin)*rho.^opt.penal;    % Reference stiffness factors
sK  = reshape(SYS.KE(:)*E',64*mesh.nelem,1);
K   = sparse(SYS.iK,SYS.jK,sK);                % Reference global stiffness
K   = 0.5*(K+K');
freeDofs = fem.bc.freeDofs;
U   = zeros(mesh.ndofs,1);
U(freeDofs) = K(freeDofs,freeDofs) \ ...
    fem.bc.forceVector(freeDofs);               % Reference displacement field
C0  = U(freeDofs)'*fem.bc.forceVector(freeDofs);% Reference compliance
[xRC,yRC] = Calculate_RC0(fem,K,U);
RC        = [xRC,yRC];                          % Reference reaction center
M0   = sum(rho.*mesh.area);                     % Initial material measure
xCG0 = sum(rho.*mesh.area.*mesh.info.Centroid(:,1))/M0;
yCG0 = sum(rho.*mesh.area.*mesh.info.Centroid(:,2))/M0;
D02  = max((xCG0-xRC)^2+(yCG0-yRC)^2,opt.tol);  % Reference squared distance
end

%----------------------------- PASSIVE PHASES
function str = enforcePassiveDensity(str,fem)
if isfield(fem.bc,'passive_void')
    str(fem.bc.passive_void) = 0;              % Prescribed void elements
end
if isfield(fem.bc,'passive_solid')
    str(fem.bc.passive_solid) = 1;             % Prescribed solid elements
end
end

%----------------------------- HISTORY ARRAYS
function h = allocateHistory(n,enabled)
if enabled
    h = struct( ...
        'iter',zeros(n,1), ...
        'J',zeros(n,1), ...
        'C',zeros(n,1), ...
        'D2',zeros(n,1), ...
        'volume',zeros(n,1));
else
    h = struct('iter',[],'J',[],'C',[],'D2',[],'volume',[]);
end
end
%-------------------------------------------------------------------------%
