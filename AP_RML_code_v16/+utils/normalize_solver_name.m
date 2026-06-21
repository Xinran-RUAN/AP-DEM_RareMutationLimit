function solver = normalize_solver_name(solverName)
%NORMALIZE_SOLVER_NAME Canonical short solver labels used in filenames.
if nargin < 1 || isempty(solverName)
    solverName = 'wkb';
end
solver = lower(strtrim(char(solverName)));
switch solver
    case {'wkb','wkb1d','wkb-1d','timeap','time-ap','projective','implicit-coupled'}
        solver = 'wkb';
    case {'direct','dir','dir-1d','direct1d','direct-1d'}
        solver = 'direct';
    case {'compare','comparison','paired','test3'}
        solver = 'compare';
    otherwise
        solver = regexprep(solver, '[^a-z0-9]+', '_');
        solver = regexprep(solver, '^_+|_+$', '');
        if isempty(solver)
            solver = 'run';
        end
end
end
