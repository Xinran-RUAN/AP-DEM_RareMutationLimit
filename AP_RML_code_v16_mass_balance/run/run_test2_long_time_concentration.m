% RUN_TEST2_LONG_TIME_CONCENTRATION
% -------------------------------------------------------------------------
% 算例 2：长时间 trait concentration。默认只跑 WKB。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

cfg.caseName = 'test2_long_time_concentration';
cfg.solvers = {'wkb'};

cfg.profile = 'baseline';
cfg.eps = 1e-2;
cfg.Nx = 50;
cfg.Ntheta_wkb = 64;
cfg.dt = 1e-3;
cfg.T = 10;
cfg.tol = 1e-8;

cfg.u0_label = 'u_flat';
cfg.W0_label = 'W_flat';
cfg.u0_fun = @(theta) 0.*theta;
cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta));

cfg.timeIntegrator = 'frozen';
cfg.phaseHamiltonian = 'weno5';
cfg.amplitudeVariant = 'split';
cfg.rhoReconstruction = 'laplace-hybrid';

% 质量平衡修正：每步用 eps*dM/dt = int rho*(K-rho) dx 修正整体质量。
cfg.massCorrectionDuringSolve = true;       % {true,false}
cfg.massCorrectionFormula = 'trapezoid';    % {'trapezoid','backward-euler'}
cfg.reactionDiscretization = 'implicit';

cfg.snapshotTimes = unique([0, 1:1:cfg.T]);
cfg.historyEveryTime = 0.2;
cfg.saveFinalResult = true;
cfg.clearOutput = true;
cfg.verbose = true;
cfg.livePlot = false;

%% ========================= 运行 =========================

runInfo = run_case_1d(cfg); %#ok<NASGU>
