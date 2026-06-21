% RUN_TIMEAP_COMPARE_1D
% -------------------------------------------------------------------------
% 时间 AP 方法对比。唯一批量变量是 methodList；本脚本逐个调用 run_case_1d。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

methodList = {'frozen', 'projective', 'implicit-coupled'};

base.caseNamePrefix = 'timeap_compare';
base.solvers = {'wkb'};

base.profile = 'baseline';
base.eps = 1e-3;
base.Nx = 40;
base.Ntheta_wkb = 32;
base.dt = 5e-3;
base.T = 1;
base.tol = 1e-8;

base.u0_label = 'u_peak_070';
base.W0_label = 'W_flat';
base.u0_fun = @(theta) -0.05*(1 - cos(2*pi*(theta - 0.70)));
base.W0_fun = @(x, theta) ones(numel(x), numel(theta));

base.phaseHamiltonian = 'weno5';
base.amplitudeVariant = 'wkb-if-split';
base.rhoReconstruction = 'laplace-hybrid';

base.massCorrectionDuringSolve = true;       % {true,false}
base.massCorrectionFormula = 'trapezoid';    % {'trapezoid','backward-euler'}
base.reactionDiscretization = 'logistic-exact';
base.hamiltonianGauge = 'min';
base.ifMode = 'h-only';
base.ifGaugeStrategy = 'phase-peak';
base.enableStepRetry = true;

base.snapshotTimes = [0, 0.5, base.T];
base.historyEveryTime = 0.1;
base.clearOutput = true;
base.verbose = true;
base.livePlot = false;

%% ========================= 运行 =========================

allRunInfo = cell(numel(methodList), 1);
for im = 1:numel(methodList)
    cfg = base;
    cfg.caseName = sprintf('%s_%s', base.caseNamePrefix, methodList{im});
    cfg.timeIntegrator = methodList{im};
    cfg.dtMax = base.dt;

    fprintf('\n========== time AP method: %s ==========%s', methodList{im}, newline);
    allRunInfo{im} = run_case_1d(cfg);
end
