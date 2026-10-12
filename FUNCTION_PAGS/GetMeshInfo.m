
%----------------------------- GetMeshInfo.m -----------------------------%
% Compute element centroids and areas for a structured Q4 mesh.
%-------------------------------------------------------------------------%

function MeshInfo = GetMeshInfo(fem)
nelx = fem.domain.nelx;
nely = fem.domain.nely;
lx = fem.domain.lx;
ly = fem.domain.ly;
x = linspace(0,lx,nelx+1);          % Nodal x-coordinates
y = linspace(0,ly,nely+1);          % Nodal y-coordinates
xc = (x(1:end-1)+x(2:end))/2;       % Element-center x-coordinates
yc = (y(1:end-1)+y(2:end))/2;       % Element-center y-coordinates
[Xc,Yc] = meshgrid(xc,yc);          % Structured element centroids
MeshInfo.Centroid = [Xc(:),Yc(:)];  % Centroid coordinates [nelem x 2]
MeshInfo.Area = (lx/nelx)*(ly/nely)* ...
                ones(nelx*nely,1);  % Uniform element areas [nelem x 1]
%-------------------------------------------------------------------------%
