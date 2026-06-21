% RUN_TEST4_COARSE_THETA_AP
% -------------------------------------------------------------------------
% 算例 4：粗 theta 网格上的 WKB AP 试验。唯一批量变量是 Ntheta。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

cfg.caseName = 'test4_coarse_theta_AP';
cfg.solvers = {'wkb'};

cfg.profile = 'baseline';
cfg.eps = 1e-3;
cfg.Nx = 40;
cfg.Ntheta_wkb = [8, 16, 32];    % 唯一批量变量
cfg.dt = 2e-3;
cfg.T = 2;
cfg.tol = 1e-8;

cfg.u0_label = 'u_peak_070';
cfg.W0_label = 'W_flat';
cfg.u0_fun = @(theta) -0.05*(1 - cos(2*pi*(theta - 0.70)));
cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta));

cfg.timeIntegrator = 'frozen';
cfg.phaseHamiltonian = 'weno5';
cfg.amplitudeVariant = 'wkb-if-split';
cfg.rhoReconstruction = 'laplace-hybrid';

% 质量平衡修正：每步用 eps*dM/dt = int rho*(K-rho) dx 修正整体质量。
cfg.massCorrectionDuringSolve = true;       % {true,false}
cfg.massCorrectionFormula = 'trapezoid';    % {'trapezoid','backward-euler'}
cfg.reactionDiscretization = 'logistic-exact';
cfg.hamiltonianGauge = 'min';
cfg.ifMode = 'h-only';
cfg.ifGaugeStrategy = 'phase-peak';
cfg.enableStepRetry = true;

cfg.snapshotTimes = [0, 0.5, 1, cfg.T];
cfg.historyEveryTime = 0.2;
cfg.clearOutput = true;
cfg.verbose = true;
cfg.livePlot = false;

%% ========================= 运行 =========================

runInfo = run_case_1d(cfg); %#ok<NASGU>
