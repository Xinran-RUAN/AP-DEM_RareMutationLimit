function progress_finish(prog, step, t, residual, info)
%PROGRESS_FINISH 打印求解结束信息。
%   与 progress_update 分开，是为了无论是否刚打印过进度，都能在最后
%   明确看到总耗时、最终步数和最终诊断量。
if nargin < 5 || isempty(info)
    info = struct();
end
if ~isfield(prog, 'enabled') || ~prog.enabled
    return;
end
wall = toc(prog.wall0);
msg = sprintf('[%s] 完成: step=%d', prog.label, step);
if isfield(prog, 'maxSteps') && isfinite(prog.maxSteps)
    msg = [msg, sprintf('/%d', prog.maxSteps)]; %#ok<AGROW>
end
msg = [msg, sprintf('  t=%s', fmt(t))]; %#ok<AGROW>
if isfield(prog, 'T') && isfinite(prog.T)
    msg = [msg, sprintf('/%s', fmt(prog.T))]; %#ok<AGROW>
end
msg = [msg, sprintf('  res=%s  elapsed=%s', fmt(residual), utils.format_seconds(wall))]; %#ok<AGROW>
msg = append_num(msg, info, 'theta_wkb', 'thetaWKB');
msg = append_num(msg, info, 'theta_direct', 'thetaDIR');
msg = append_num(msg, info, 'err_wkb_to_m', 'errWKB');
msg = append_num(msg, info, 'err_direct_to_m', 'errDIR');
msg = append_num(msg, info, 'gammaR', 'gammaR');
msg = append_num(msg, info, 'rho_max', 'rhoMax');
fprintf('%s\n', msg);
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
