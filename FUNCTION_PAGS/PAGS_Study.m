
%----------------------------- PAGS_Study.m ------------------------------%
% Execute and organize the numerical studies for PAGS-LSRD optimization.
% Input : fem          - Finite element model and boundary conditions
%         opt          - Optimization and material parameters
%         study        - Study configurations and output directory
% Output: studyResults - Baseline, convergence, and parametric results
%-------------------------------------------------------------------------%

function studyResults = PAGS_Study(fem,opt,study)
%--------- Validate the finite element model and execution mode
if ~isfield(fem,'bc')
    error('PAGS_Study:MissingBC','Call ApplyBC first.');
end
if ~isfield(study,'runMode')
    error('PAGS_Study:MissingMode','Missing runMode.');
end
mode = lower(char(study.runMode));             % Selected execution mode
if ~ismember(mode,{'all','history','sweep'})
    error('PAGS_Study:InvalidMode','Use all, history, or sweep.');
end
if ~exist(study.outputFolder,'dir')
    mkdir(study.outputFolder);                 % Create output directory
end

%--------- Initialize mesh 
mesh = struct('fem',fem,'info',GetMeshInfo(fem));
mesh.nelx = fem.domain.nelx;
mesh.nely = fem.domain.nely;
mesh.nelem = mesh.nelx*mesh.nely;            % Number of Q4 elements
mesh.ndofs = 2*(mesh.nelx+1)*(mesh.nely+1);  % Number of displacement DOFs
mesh.area = mesh.info.Area(:);               % Element area vector

[mesh.X,mesh.Y] = meshgrid(linspace(0,fem.domain.lx,mesh.nelx+1), ...
    linspace(0,fem.domain.ly,mesh.nely+1));     % Nodal coordinates for plotting

studyResults = struct('history',[],'sweep',[]); % Initialize result storage

%----------------------------- CONVERGENCE STUDY
% Reproduce the reference and orientation-dependent results (Figs. 3-6).
if ismember(mode,{'all','history'})
    cfg  = study.history;
    optH = opt;
    optH.alpha = cfg.alpha;                     % Level-set relaxation
    optH.delta = cfg.delta;                     % Heaviside half-width

    baseCase = makeCase(0,cfg.referenceW,optH.penal,true,true);
    baseline = PAGS_LSRD(mesh,optH,baseCase);   % Compliance-only reference
    previousReferencePenal = optH.penal;
    for ia = 1:numel(cfg.theta)
        caseDef = makeCase(cfg.theta(ia),cfg.w,previousReferencePenal, ...
                                                               true,false);
        one = PAGS_LSRD(mesh,optH,caseDef); % Orientation-dependent case
        if ia == 1
            angleCases = one;
        else
            angleCases(ia) = one;
        end
        previousReferencePenal = one.finalPenal;
    end

    history = struct('baseline',baseline,'angles',angleCases);
    studyResults.history = history;            % Store convergence results
    PAGS_Output('history',history,mesh,study.outputFolder);
end

%----------------------------- PARAMETRIC STUDY
if ismember(mode,{'all','sweep'})
    cfg = study.sweep;
    optS = opt;
    optS.alpha = cfg.alpha;                    % Level-set relaxation
    optS.delta = cfg.delta;                    % Heaviside half-width

    for ia = 1:numel(cfg.theta)
        for iw = 1:numel(cfg.w)

            caseDef = makeCase(cfg.theta(ia),cfg.w(iw), ...
                               1.0,false,false);
            one = PAGS_LSRD(mesh,optS,caseDef); % Independent parameter case

            if ia == 1 && iw == 1
                cases = repmat(one,numel(cfg.theta),numel(cfg.w));
            end
            cases(ia,iw) = one;                 % Store case results
        end
    end
    studyResults.sweep = cases;
    PAGS_Output('sweep',cases,mesh,study.outputFolder);
end
fprintf('Published numerical outputs: %s\n',study.outputFolder);
end

%----------------------------- CASE DEFINITION
function def = makeCase(thetaDeg,w,referencePenal,recordHistory,recordSnapshots)
def = struct( ...
    'thetaDeg',thetaDeg, ...
    'w',w, ...
    'referencePenal',referencePenal, ...
    'recordHistory',recordHistory, ...
    'recordSnapshots',recordSnapshots);
end
%-------------------------------------------------------------------------%
