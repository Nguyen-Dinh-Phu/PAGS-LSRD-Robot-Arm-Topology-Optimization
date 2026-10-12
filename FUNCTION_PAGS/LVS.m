
%-------------------------------- LVS.m ----------------------------------%
% Update the level-set field using reaction-diffusion evolution with
% volume continuation and an exponential Lagrange multiplier.
%-------------------------------------------------------------------------%

function [phi,str] = LVS(phi,TDN,vol,volInit,iterNum,nvol,Vmax, ...
                         betaVol,dVol,NN,dt,T,freedofsphi,nelx,nely)
%--------- Prescribed volume continuation
if iterNum <= nvol
targetVol     = Vmax + (volInit-Vmax)*(1-iterNum/nvol);
else
    targetVol = Vmax;
end
targetVol     = max(targetVol, 1e-10);              % Positive target volume

%--------- Nodal sensitivity normalization
g = reshape(TDN,(nely+1)*(nelx+1),1);          % Nodal driving field
meanSens = max(mean(abs(g)),1e-12);            % Sensitivity scaling
Cscale   = 1/meanSens;                         % Normalization factor

%--------- Exponential volume constraint multiplier
volumeError = (vol-targetVol)/targetVol;       % Relative volume deviation
lambda      = meanSens*exp(betaVol*(volumeError+dVol));

%--------- Reaction-diffusion level-set update
reaction = Cscale*(g-lambda);                  % Sensitivity-driven reaction
rhs      = NN*(reaction+phi/dt);               % Discrete evolution RHS

phiNew = phi;
phiNew(freedofsphi) = T(freedofsphi,freedofsphi) \ rhs(freedofsphi);
phi    = min(1,max(-1,phiNew));                % Bound level-set values

%--------- Reconstruct element densities
str    = PhiToDensity(phi,nelx,nely,0.10);     % Regularized Heaviside mapping
%-------------------------------------------------------------------------%
