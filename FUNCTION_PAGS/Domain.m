
%------------------------------- Domain.m --------------------------------%
% Generate a structured four-node quadrilateral (Q4) finite element mesh.
% Input : lx, ly     - Domain dimensions
%         nelx, nely - Number of elements along x and y
% Output: fem        - Domain geometry and mesh connectivity
%-------------------------------------------------------------------------%

function fem = Domain(lx,ly,nelx,nely)
fem.domain = struct('lx',lx,'ly',ly,'nelx',nelx,'nely',nely);
[X,Y] = meshgrid(linspace(0,lx,nelx+1),linspace(0,ly,nely+1));
fem.mesh.node = [X(:),Y(:)];            % Nodal coordinates
[i,j] = meshgrid(1:nelx,1:nely);
n1 = (i(:)-1)*(nely+1) + j(:);          % Lower-left nodes
n3 = n1 + nely + 1;                     % Lower-right nodes
fem.mesh.elem = [n1,n3,n3+1,n1+1];      % Counterclockwise Q4 connectivity
%-------------------------------------------------------------------------%
