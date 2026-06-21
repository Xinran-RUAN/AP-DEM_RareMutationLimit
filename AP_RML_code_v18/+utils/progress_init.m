function prog = progress_init(label, par)
%PROGRESS_INIT 初始化主程序/求解器的进度输出状态。
%   这个函数只负责记录开始时间、预计总步数和输出频率。真正的
%   进度打印由 utils.progress_update 完成。把这部分做成工具函数，
%   可以保证 WKB、direct 和各个 run 脚本的输出格式一致。
if nargin < 1 || isempty(label)
    label = 'RUN';
end
if nargin < 2 || isempty(par)
    par = struct();
end

prog = struct();
prog.label = label;
prog.enabled = true;
if isfield(par, 'verbose') && ~isempty(par.verbose)
    prog.enabled = logical(par.verbose);
end
prog.wall0 = tic;
prog.lastWall = -Inf;
prog.lastStep = -Inf;
prog.nextPercent = 0;

prog.T = get_num(par, 'T', NaN);
prog.dt = get_num(par, 'dt', NaN);
prog.tol = get_num(par, 'tol', -Inf);
prog.stopByResidual = get_bool(par, 'stopByResidual', false);
prog.progressPrintWhenTolReached = get_bool(par, 'progressPrintWhenTolReached', false);
maxStepsInput = get_num(par, 'maxSteps', Inf);
if isfinite(prog.T) && isfinite(prog.dt) && prog.dt > 0
    stepsByTime = ceil(prog.T / prog.dt);
else
    stepsByTime = maxStepsInput;
end
if isfinite(maxStepsInput)
    prog.maxSteps = min(stepsByTime, maxStepsInput);
else
    prog.maxSteps = stepsByTime;
end
if ~(isfinite(prog.maxSteps) && prog.maxSteps > 0)
    prog.maxSteps = NaN;
end

prog.everyPercent = get_num(par, 'progressEveryPercent', 5);
if ~(isfinite(prog.everyPercent) && prog.everyPercent > 0)
    prog.everyPercent = 5;
end
prog.everySeconds = get_num(par, 'progressEverySeconds', 10);
if ~(isfinite(prog.everySeconds) && prog.everySeconds >= 0)
    prog.everySeconds = 10;
end
prog.everySteps = get_num(par, 'progressEverySteps', NaN);
if ~(isfinite(prog.everySteps) && prog.everySteps > 0)
    if isfinite(prog.maxSteps)
        prog.everySteps = max(1, ceil(prog.maxSteps * prog.everyPercent / 100));
    else
        prog.everySteps = 100;
    end
end

if prog.enabled
    fprintf('\n[%s] 开始计算\n', prog.label);
    fprintf('  参数: eps=%s, Nx=%s, Ntheta=%s, T=%s, dt=%s, tol=%s\n', ...
        fmt(get_num(par,'eps',NaN)), fmt(get_num(par,'Nx',NaN)), ...
        fmt(get_num(par,'Ntheta',NaN)), fmt(prog.T), fmt(prog.dt), fmt(prog.tol));
    if isfinite(prog.maxSteps)
        fprintf('  预计最多步数: %d；进度输出: 每 %.1f%% 或 %.1fs 或 %d 步\n', ...
            prog.maxSteps, prog.everyPercent, prog.everySeconds, prog.everySteps);
    else
        fprintf('  预计最多步数: 未知；进度输出: 每 %.1fs 或 %d 步\n', ...
            prog.everySeconds, prog.everySteps);
    end
    if isfield(par, 'outdir') && ~isempty(par.outdir)
        fprintf('  输出目录: %s\n', char(par.outdir));
    end
end
end

function val = get_num(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    tmp = s.(name);
    if isnumeric(tmp) && isscalar(tmp)
        val = double(tmp);
    end
end
end


function val = get_bool(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = logical(s.(name));
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
