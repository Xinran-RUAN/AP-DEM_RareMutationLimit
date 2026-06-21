function prog = progress_update(prog, step, t, residual, info)
%PROGRESS_UPDATE 按固定步数/百分比/墙钟时间打印进度。
%   info 可以包含 theta_wkb、theta_direct、err_to_m、gammaR、rho_max、
%   sigma_theta 等字段；没有的字段会自动跳过。
if nargin < 5 || isempty(info)
    info = struct();
end
if nargin < 4 || isempty(residual)
    residual = NaN;
end
if nargin < 3 || isempty(t)
    t = NaN;
end
if ~isfield(prog, 'enabled') || ~prog.enabled
    return;
end

wall = toc(prog.wall0);
force = false;
if isfield(info, 'force') && ~isempty(info.force) && isscalar(info.force)
    force = logical(info.force);
end

pct = NaN;
if isfield(prog, 'maxSteps') && isfinite(prog.maxSteps) && prog.maxSteps > 0
    pct = min(100, 100 * step / prog.maxSteps);
end

hitPercent = isfinite(pct) && pct >= prog.nextPercent;
hitStep = step <= 1 || step - prog.lastStep >= prog.everySteps;
hitWall = wall - prog.lastWall >= prog.everySeconds;
hitEnd = false;
if isfield(prog, 'maxSteps') && isfinite(prog.maxSteps)
    hitEnd = step >= prog.maxSteps;
end
if isfield(prog, 'T') && isfinite(prog.T)
    hitEnd = hitEnd || t >= prog.T - 1e-14;
end
hitTol = isfield(prog, 'tol') && isfinite(prog.tol) && residual <= prog.tol;
% 默认只有 stopByResidual=true 或显式 progressPrintWhenTolReached=true 时，
% 才因为 residual<=tol 强制打印。否则收敛后每一步都会命中 hitTol。
printTol = hitTol && ( ...
    (isfield(prog,'stopByResidual') && logical(prog.stopByResidual)) || ...
    (isfield(prog,'progressPrintWhenTolReached') && logical(prog.progressPrintWhenTolReached)) );

if ~(force || hitPercent || hitStep || hitWall || hitEnd || printTol)
    return;
end

rate = step / max(wall, 1e-12);
eta = NaN;
if isfield(prog, 'maxSteps') && isfinite(prog.maxSteps) && rate > 0
    eta = max(prog.maxSteps - step, 0) / rate;
end

if isfinite(pct)
    msg = sprintf('[%s] %6d/%6d (%5.1f%%)', prog.label, step, prog.maxSteps, pct);
else
    msg = sprintf('[%s] step=%6d', prog.label, step);
end
msg = [msg, sprintf('  t=%s', fmt(t))]; %#ok<AGROW>
if isfield(prog, 'T') && isfinite(prog.T)
    msg = [msg, sprintf('/%s', fmt(prog.T))]; %#ok<AGROW>
end
msg = [msg, sprintf('  res=%s  elapsed=%s  ETA=%s', ...
    fmt(residual), utils.format_seconds(wall), utils.format_seconds(eta))]; %#ok<AGROW>

msg = append_num(msg, info, 'theta', 'theta');
msg = append_num(msg, info, 'theta_wkb', 'thetaWKB');
msg = append_num(msg, info, 'theta_direct', 'thetaDIR');
msg = append_num(msg, info, 'err_to_m', 'err');
msg = append_num(msg, info, 'err_wkb_to_m', 'errWKB');
msg = append_num(msg, info, 'err_direct_to_m', 'errDIR');
msg = append_num(msg, info, 'sigma_theta', 'sigma');
msg = append_num(msg, info, 'gammaR', 'gammaR');
msg = append_num(msg, info, 'rho_max', 'rhoMax');
msg = append_num(msg, info, 'Hmin', 'Hmin');
msg = append_num(msg, info, 'dt', 'dt');
msg = append_num(msg, info, 'res_n', 'resN');
msg = append_num(msg, info, 'res_u', 'resU');
msg = append_num(msg, info, 'res_W_legacy', 'resWold');
msg = append_num(msg, info, 'res_phase_steady', 'resPhase');
msg = append_num(msg, info, 'res_density_steady', 'resDens');
msg = append_num(msg, info, 'res_amp_w_steady', 'resAmpW');
msg = append_num(msg, info, 'lambdaR', 'lambdaR');
msg = append_num(msg, info, 'lambdaH', 'lambdaH');
msg = append_num(msg, info, 'phaseCFL', 'phaseCFL');
msg = append_num(msg, info, 'ifLogRangeEstimate', 'ifLogR');
msg = append_num(msg, info, 'jumpEps', 'jump/eps');
msg = append_num(msg, info, 'peakGrid', 'peakGrid');
msg = append_num(msg, info, 'logQMax', 'logQmax');
msg = append_num(msg, info, 'ifGaugeShift', 'gShift');
msg = append_num(msg, info, 'numRetries', 'retry');
msg = append_num(msg, info, 'seedLogRange', 'seedRange');
msg = append_num(msg, info, 'seedLogMax', 'seedMax');
msg = append_num(msg, info, 'reactionScaleMax', 'reactScale');
msg = append_num(msg, info, 'rhoLogMax', 'logRho');
msg = append_num(msg, info, 'logWRangeAfter', 'logWrange');
msg = append_num(msg, info, 'thetaRephaseARange', 'rephaseA');
msg = append_num(msg, info, 'implicitIter', 'impIt');
msg = append_num(msg, info, 'implicitRelChange', 'impRel');
msg = append_num(msg, info, 'projectiveMicroStepsUsed', 'microN');
msg = append_num(msg, info, 'projectiveMicroDt', 'microDt');
if isfield(info, 'note') && ~isempty(info.note)
    msg = [msg, '  ', char(info.note)]; %#ok<AGROW>
end
fprintf('%s\n', msg);

prog.lastWall = wall;
prog.lastStep = step;
if isfinite(pct)
    while pct >= prog.nextPercent
        prog.nextPercent = prog.nextPercent + prog.everyPercent;
    end
end
end

function msg = append_num(msg, s, field, label)
if isstruct(s) && isfield(s, field) && ~isempty(s.(field)) && isnumeric(s.(field)) && isscalar(s.(field))
    x = double(s.(field));
    if isfinite(x)
        msg = [msg, sprintf('  %s=%s', label, fmt(x))]; %#ok<AGROW>
    end
end
end

function str = fmt(x)
if isempty(x) || ~isfinite(x)
    str = '--';
elseif abs(x) >= 1e4 || (abs(x) > 0 && abs(x) < 1e-3)
    str = sprintf('%.3e', x);
else
    str = sprintf('%.6g', x);
end
end
