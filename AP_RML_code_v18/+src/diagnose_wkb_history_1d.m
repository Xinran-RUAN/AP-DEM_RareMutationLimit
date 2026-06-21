function report = diagnose_wkb_history_1d(history, doPrint)
%DIAGNOSE_WKB_HISTORY_1D  中文诊断 WKB 小 epsilon 计算是否病态。
%
% 用法：
%   report = src.diagnose_wkb_history_1d(wkb.history, true);
%
% 重点查看：
%   res_phase_steady      相位 HJ 稳态残差；
%   res_amplitude_steady  W 振幅稳态残差；
%   qLogAbsMax            IF 指数大小；
%   seedLogRange/Max      qW 的 log 尺度，过大时会拒步缩小 dt；
%   uPeakStar             phase-peak gauge 后应接近 0；
%   peakWidthOverGrid     密度峰宽/trait 网格宽度，小于 1 表示密度峰形欠解析；
%   numRetries            是否发生拒步重试。

if nargin < 2, doPrint = true; end
report = struct();
if isempty(history) || ~isstruct(history)
    report.message = 'history 为空，无法诊断。';
    if doPrint, fprintf('%s\n', report.message); end
    return;
end

fields = fieldnames(history);
getlast = @(name) local_last(history, name);
report.stepsRecorded = numel(getfield(history, fields{1})); %#ok<GFLD>
report.tLast = getlast('t');
report.residualLast = getlast('residual');
report.resPhaseLast = getlast('res_phase_steady');
report.resAmpLast = getlast('res_amplitude_steady');
report.resDensLast = getlast('res_density_steady');
report.dtLast = getlast('dt');
report.phaseCFLLast = getlast('phaseCFL');
report.qLogAbsMaxLast = getlast('qLogAbsMax');
report.seedLogRangeLast = getlast('seedLogRange');
report.seedLogMaxLast = getlast('seedLogMax');
report.numRetriesMax = local_max(history, 'numRetries');
report.peakWidthOverGridLast = getlast('peakWidthOverGrid');
report.uPeakStarLast = getlast('uPeakStar');
report.rhoMaxLast = getlast('rho_max');
report.rhoLogMaxLast = getlast('rhoReconstructLogMax');
report.logWRangeLast = getlast('logWRange');
report.thetaRephaseARangeLast = getlast('thetaRephaseARange');
report.implicitIterLast = getlast('implicitIter');
report.implicitRelChangeLast = getlast('implicitRelChange');
report.implicitConvergedLast = getlast('implicitConverged');
report.projectiveMicroStepsLast = getlast('projectiveMicroStepsUsed');
report.projectiveMicroDtLast = getlast('projectiveMicroDt');

messages = {};
if isfinite(report.peakWidthOverGridLast) && report.peakWidthOverGridLast < 1
    messages{end+1} = '密度峰宽小于一个 trait 网格：thetaDIR/密度峰形不可靠，应主要看 WKB 相位峰位置。'; %#ok<AGROW>
end
if isfinite(report.qLogAbsMaxLast) && report.qLogAbsMaxLast > 80
    messages{end+1} = 'IF 指数 q 的绝对值过大：建议降低 ifLogRangeMax 或查看 H 的跨度。'; %#ok<AGROW>
end
if isfinite(report.seedLogRangeLast) && report.seedLogRangeLast > 160
    messages{end+1} = 'qW 的 log 跨度过大：这是 W 振幅变量中的主要病态来源。若 v9 的 theta 再分解开启后仍很大，需要检查 H 跨度或降低 ifLogRangeMax。'; %#ok<AGROW>
end
if isfinite(report.logWRangeLast) && report.logWRangeLast > 100
    messages{end+1} = 'W 的 log 范围仍然较大：theta 方向再分解没有充分吸收指数尺度，可尝试增加 thetaRephaseSmoothPasses 或检查 W 是否出现局部污染。'; %#ok<AGROW>
end
if isfinite(report.uPeakStarLast) && abs(report.uPeakStarLast) > 1e-8
    messages{end+1} = 'phase-peak gauge 后 uPeakStar 没有接近 0，需检查 gauge 或 theta 再分解。'; %#ok<AGROW>
end
if isfinite(report.implicitRelChangeLast) && report.implicitRelChangeLast > 1e-3
    messages{end+1} = sprintf('全隐式 Picard 最后相对变化 %.3g，说明 H^{n+1} 耦合还没有充分收敛；可减小 implicitRelax 或增加 implicitMaxIter。', report.implicitRelChangeLast); %#ok<AGROW>
end
if isfinite(report.projectiveMicroDtLast) && report.projectiveMicroDtLast > 0
    messages{end+1} = sprintf('projective 模式最后 microDt=%.3e，micro steps=%.0f；若结果不稳，可增加 micro steps 或减小 projectiveMicroDtFactor。', report.projectiveMicroDtLast, report.projectiveMicroStepsLast); %#ok<AGROW>
end
if isfinite(report.numRetriesMax) && report.numRetriesMax > 0
    messages{end+1} = sprintf('曾发生拒步重试，最大重试次数为 %.0f；说明原 dt 对某些瞬间过大。', report.numRetriesMax); %#ok<AGROW>
end
if isfinite(report.resPhaseLast) && isfinite(report.resAmpLast)
    if report.resPhaseLast > 10*report.resAmpLast
        messages{end+1} = '当前 residual 主要来自相位方程；可检查 phaseCFL 或 H 计算。'; %#ok<AGROW>
    elseif report.resAmpLast > 10*report.resPhaseLast
        messages{end+1} = '当前 residual 主要来自 W 振幅方程；可检查 qW log 尺度、reactionScale 和 W 矩阵。'; %#ok<AGROW>
    end
end
if isempty(messages)
    messages{1} = '未发现明显单一病态来源；请同时查看各分量曲线。';
end
report.messages = messages;

if doPrint
    fprintf('\n========== WKB 小 epsilon 中文诊断 ==========' );
    fprintf('\n记录点数: %d, 最后时间 t=%.6g, 最后 residual=%.6g\n', report.stepsRecorded, report.tLast, report.residualLast);
    fprintf('resPhase=%.6g, resAmp=%.6g, resDens=%.6g, dt=%.3e, phaseCFL=%.3g\n', ...
        report.resPhaseLast, report.resAmpLast, report.resDensLast, report.dtLast, report.phaseCFLLast);
    fprintf('qLogAbsMax=%.3g, seedLogMax=%.3g, seedLogRange=%.3g, maxRetries=%.0f\n', ...
        report.qLogAbsMaxLast, report.seedLogMaxLast, report.seedLogRangeLast, report.numRetriesMax);
    fprintf('peakWidth/grid=%.3g, uPeakStar=%.3e, rhoMax=%.3g, logRhoMax=%.3g\n', ...
        report.peakWidthOverGridLast, report.uPeakStarLast, report.rhoMaxLast, report.rhoLogMaxLast);
    fprintf('logWRange=%.3g, thetaRephaseARange=%.3g\n', report.logWRangeLast, report.thetaRephaseARangeLast);
    fprintf('implicitIter=%.0f, implicitRel=%.3g, implicitConv=%.0f, microN=%.0f, microDt=%.3e\n', ...
        report.implicitIterLast, report.implicitRelChangeLast, report.implicitConvergedLast, ...
        report.projectiveMicroStepsLast, report.projectiveMicroDtLast);
    for i = 1:numel(messages)
        fprintf('诊断 %d: %s\n', i, messages{i});
    end
    fprintf('============================================\n');
end
end

function v = local_last(h, name)
if isfield(h,name) && ~isempty(h.(name))
    a = h.(name);
    v = a(end);
else
    v = NaN;
end
end

function v = local_max(h, name)
if isfield(h,name) && ~isempty(h.(name))
    a = h.(name);
    v = max(a(isfinite(a)));
    if isempty(v), v = NaN; end
else
    v = NaN;
end
end
