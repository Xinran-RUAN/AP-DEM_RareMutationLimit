function par = apply_patch_preset_1d(par, preset)
%APPLY_PATCH_PRESET_1D 为合并版一维代码应用一组“显式可选”的补丁预设。
%
%  用法示例：
%      par = model.default_params_1d('baseline');      % 旧版默认
%      par = model.apply_patch_preset_1d(par,'v13-small-eps');
%
%  preset 可取：
%    'legacy'          : 强制恢复旧版默认，适合复现实验；
%    'v13-small-eps'   : 打开 v13 小 epsilon 稳健 WKB-IF 设置；
%    'v13-projective'  : 在 small-eps 基础上启用 projective 时间 AP 实验；
%    'v13-implicit'    : 在 small-eps 基础上启用 implicit-coupled 时间 AP 实验。
%
%  注意：本函数只改 par，不直接运行求解器；所有补丁都可在调用后继续手动覆盖。

if nargin < 2 || isempty(preset)
    preset = 'legacy';
end
preset = lower(char(preset));

switch preset
    case {'legacy','old','paper'}
        par.patchPreset = 'legacy';
        par.initialPhaseMode = 'flat';
        par.initialAmplitudeMode = 'flat-rho';
        par.initialGaugeMode = 'legacy-max';
        par.prepareInitialMass = false;
        par.timeIntegrator = 'frozen';
        par.phaseHamiltonian = 'weno5';
        par.amplitudeVariant = 'split';
        par.reactionDiscretization = 'implicit';
        par.rhoReconstruction = 'direct-log';
        par.residualMode = 'legacy';
        par.useSubgridPhasePeak = false;
        par.smallEpsAutoPatches = false;
        par.enableV13AutoPatches = false;
        par.enableSmallEpsMonitor = false;
        par.hamiltonianGauge = 'none';
        par.ifMode = 'none';
        par.ifGaugeStrategy = 'none';
        par.gaugeMode = 'max';
        par.postGaugeMode = 'max';
        par.thetaRephase = false;
        par.adaptiveTimeStep = false;
        par.adaptiveIfLogRange = false;
        par.enableStepRetry = false;
        par.rhoLogRetryMax = Inf;

    case {'v13-small-eps','small-eps','v13','new','wkb-if'}
        par.patchPreset = 'v13-small-eps';
        % 初值仍可保留 flat；这里不强迫改成 peaked，避免和旧实验物理设定混淆。
        % 若需要 v13 曾用的 peaked+well-prepared 初值，可在调用后设置：
        %   par.initialPhaseMode='peaked'; par.prepareInitialMass=true;
        par.initialGaugeMode = 'phase-peak';
        par.prepareInitialMass = true;
        par.initialRhoProfile = 'scaledK';
        par.initialRhoScale = 0.5;
        par.timeIntegrator = 'frozen';
        par.phaseHamiltonian = 'godunov';
        par.amplitudeVariant = 'wkb-if-split';
        par.reactionDiscretization = 'logistic-exact';
        par.rhoReconstruction = 'laplace-hybrid';
        par.residualMode = 'wkb-steady';
        par.useSubgridPhasePeak = true;
        par.smallEpsAutoPatches = true;
        par.enableV13AutoPatches = true;
        par.enableSmallEpsMonitor = true;
        par.hamiltonianGauge = 'min';
        par.ifMode = 'h-only';
        par.ifGaugeStrategy = 'phase-peak';
        par.gaugeMode = 'auto';
        par.postGaugeMode = 'none';
        par.thetaRephase = false;
        par.adaptiveIfLogRange = true;
        par.ifLogRangeMax = 80;
        par.enableStepRetry = true;
        par.maxStepRetries = 12;
        par.retryDtFactor = 0.5;
        par.rhoLogRetryMax = 120;

    case {'v13-projective','projective','timeap-projective'}
        par = model.apply_patch_preset_1d(par, 'v13-small-eps');
        par.patchPreset = 'v13-projective';
        par.timeIntegrator = 'projective';
        par.projectiveMicroSteps = 8;
        par.projectiveMicroDtFactor = 0.25;
        par.projectiveMode = 'u-only';
        par.projectiveCorrectorSteps = 2;
        par.adaptiveTimeStep = true;
        par.adaptiveStrategy = 'phase-cfl';
        par.dtMax = par.dt;

    case {'v13-implicit','implicit','implicit-coupled','timeap-implicit'}
        par = model.apply_patch_preset_1d(par, 'v13-small-eps');
        par.patchPreset = 'v13-implicit';
        par.timeIntegrator = 'implicit-coupled';
        par.implicitMaxIter = 100;
        par.implicitTol = 1e-7;
        par.implicitRelax = 0.6;
        par.implicitFailureTol = 5e-2;
        par.adaptiveTimeStep = true;
        par.adaptiveStrategy = 'phase-cfl';
        par.dtMax = par.dt;

    otherwise
        error('Unknown patch preset "%s".', preset);
end
end
