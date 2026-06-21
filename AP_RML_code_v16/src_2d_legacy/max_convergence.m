% Legacy 2-D convergence check with canonical snapshot paths.

close all
clear

t = 2;
eps = 1e-0;
ol = 0; %#ok<NASGU>
oll = 0;

Nth = [15, 20, 50, 100, 200, 500, 1000];
Nx = 5;
x = linspace(0, 1, Nx + 1);
x = x';
M = zeros(1, length(Nth));

wkbRoot = utils.output_dir('legacy_main_bio');
wkbBaseDir = utils.find_case_dir(wkbRoot, sprintf('eps_%s', utils.time_token(eps)), 'eps_*', true);
for kk = 1:length(Nth)
    theta = linspace(0, 1, Nth(kk) + 1);
    theta(end) = [];
    theta_f = 0:0.0001:1;

    [wkbFile, ~] = utils.find_legacy_snapshot_file(wkbBaseDir, 'legacy_wkb', [], Nth(kk), t, true);
    Swkb = load(wkbFile);
    u = Swkb.u;
    W = Swkb.W;
    if isfield(Swkb, 'theta')
        theta = Swkb.theta;
    end
    if isfield(Swkb, 'eps')
        eps = Swkb.eps;
    end

    [u_inte, W_inte] = uW_spline_theta(u, W, x, theta, theta_f);
    n_eps = W_inte .* exp(u_inte/eps);
    n_theta = n_int_x(n_eps, x, theta_f, oll);
    M(kk) = max(n_theta);
end

No = [20, 50, 100, 200, 500, 1000];
Mo = zeros(1, length(No));
directRoot = utils.output_dir('legacy_original_bio');
directBaseDir = utils.find_case_dir(directRoot, sprintf('eps_%s', utils.time_token(eps)), 'eps_*', true);
for k = 1:length(No)
    N = No(k);
    theta_o = linspace(0, 1, N + 1);
    theta_o(end) = [];

    [directFile, ~] = utils.find_legacy_snapshot_file(directBaseDir, 'legacy_direct', Nx, N, t, true);
    Sdir = load(directFile);
    ne = Sdir.ne;
    ne = [ne, ne(:, 1)];
    netheta = n_int_x(ne, x, theta_o, oll);
    Mo(k) = max(netheta);
end

figDir = utils.figure_dir(wkbBaseDir);
figure
plot(log10(Nth(1:end-1)), log10(M(1:end-1) - M(end)), 'b-*');
hold on
plot(log10(No), log10(Mo - M(end)), 'r-->');
saveas(gcf, fullfile(figDir, 'legacy_max_convergence_against_wkb.png'));

figure
plot(log10(Nth), log10(M - Mo(end)), 'b-*');
hold on
plot(log10(No(1:end-1)), log10(Mo(1:end-1) - Mo(end)), 'r-->');
saveas(gcf, fullfile(figDir, 'legacy_max_convergence_against_direct.png'));
