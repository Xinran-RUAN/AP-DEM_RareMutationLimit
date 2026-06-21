function tag = build_tag(par)
%BUILD_TAG Build a compact filename tag.
amp = get_string(par, 'amplitudeVariant', 'wkb');
rho = get_string(par, 'rhoReconstruction', 'rho');
react = get_string(par, 'reactionDiscretization', '');
resmode = get_string(par, 'residualMode', '');
adapt = '';
timePart = '';
ti = get_string(par, 'timeIntegrator', '');
if ~isempty(ti) && ~strcmpi(ti,'frozen')
    timePart = ['_' ti];
end
if isfield(par, 'adaptiveTimeStep') && ~isempty(par.adaptiveTimeStep) && logical(par.adaptiveTimeStep)
    adapt = '_adapt';
end
if isempty(react)
    reactPart = '';
else
    reactPart = ['_' react];
end
if isempty(resmode)
    resPart = '';
else
    resPart = ['_' resmode];
end
tag = sprintf('%s_eps%.3g_Nx%d_K%d_dt%.3g_%s%s_%s%s%s%s', ...
    par.profile, par.eps, par.Nx, par.Ntheta, par.dt, amp, reactPart, rho, resPart, adapt, timePart);
tag = regexprep(tag, '[^A-Za-z0-9_\.-]', '_');
end

function val = get_string(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end
