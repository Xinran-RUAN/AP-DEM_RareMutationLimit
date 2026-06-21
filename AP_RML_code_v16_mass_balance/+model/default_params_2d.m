function par = default_params_2d()
%DEFAULT_PARAMS_2D 二维 legacy FEM 实验的默认参数。
%   二维核心仍沿用旧 FEM 代码，这里只集中管理参数、路径和进度输出。
par.Ntheta = 32;
par.dt = 1e-3;
par.T = 5;
par.maxSteps = Inf;
par.tol = 1e-9;
par.eps = 1e-2;
par.h = 0.08;
par.theta_m = 0.35;
par.theta0 = 0.70;
par.profile = 'baseline2d';
par.outdir = utils.output_dir('test6_2d');
par.verbose = true;
par.progressEveryPercent = 5;
par.progressEverySeconds = 10;
par.progressEverySteps = [];
par.historyEveryTime = 0.1;
par.historyEverySteps = [];

theta_m = par.theta_m;
par.K_fun = @(x,y) 1 + 0.35*cos(2*pi*x).*cos(2*pi*y) + 0.15*cos(4*pi*x);
par.D_fun = @(theta) 0.2 + 0.4*(1 - cos(2*pi*(theta - theta_m)));
end
