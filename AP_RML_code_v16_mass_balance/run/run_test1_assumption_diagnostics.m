% RUN_TEST1_ASSUMPTION_DIAGNOSTICS
% -------------------------------------------------------------------------
% 算例 1：不同 eps 下的基本诊断。唯一批量变量是 eps。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

cfg.caseName = 'test1_assumption_diagnostics';
cfg.solvers = {'wkb'};

cfg.profile = 'baseline';
cfg.eps = [1e-1, 1e-2];       % 唯一批量变量
cfg.Nx = 40;
cfg.Ntheta_wkb = 64;
cfg.dt = 1e-3;
cfg.T = 2;
cfg.tol = 1e-8;

cfg.u0_label = 'u_flat';
cfg.W0_label = 'W_flat';
cfg.u0_fun = @(theta) 0.*theta;
cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta));

cfg.timeIntegrator = 'frozen';
cfg.phaseHamiltonian = 'weno5';
cfg.amplitudeVariant = 'split';
cfg.rhoReconstruction = 'laplace-hybrid';
cfg.reactionDiscretization = 'implicit';
cfg.residualMode = 'legacy';

% 质量平衡修正：每步用 eps*dM/dt = int rho*(K-rho) dx 修正整体质量。
cfg.massCorrectionDuringSolve = true;       % {true,false}
cfg.massCorrectionFormula = 'trapezoid';    % {'trapezoid','backward-euler'}

cfg.snapshotTimes = [0, 1, cfg.T];
cfg.historyEveryTime = 0.2;
cfg.clearOutput = true;
cfg.verbose = true;
cfg.livePlot = false;

%% ========================= 运行 =========================

runInfo = run_case_1d(cfg); %#ok<NASGU>
