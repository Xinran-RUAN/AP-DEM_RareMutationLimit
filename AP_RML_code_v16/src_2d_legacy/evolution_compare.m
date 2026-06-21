% Compare legacy 2-D WKB reconstruction with legacy direct density.
% All paths/names are canonical; no ./data or ./original_data hard coding.

close all
clear

t = 2;
eps = 1e-1;
ol = 0; %#ok<NASGU>
oll = 1;

Nth = 15;
Nx = 5;
theta = linspace(0, 1, Nth + 1);
theta(end) = [];
x = linspace(0, 1, Nx + 1);
x = x';
theta_f = 0:0.001:1;

wkbRoot = utils.output_dir('legacy_main_bio');
wkbDir = utils.find_case_dir(wkbRoot, sprintf('eps_%s', utils.time_token(eps)), 'eps_*', true);
[wkbFile, wkbTime] = utils.find_legacy_snapshot_file(wkbDir, 'legacy_wkb', [], Nth, t, true);
Swkb = load(wkbFile);
u = Swkb.u;
W = Swkb.W;
if isfield(Swkb, 'theta')
    theta = Swkb.theta;
end
if isfield(Swkb, 'eps')
    eps = Swkb.eps;
end
fprintf('[legacy compare] WKB snapshot: %s, t=%g\n', wkbFile, wkbTime);

[u_inte, W_inte] = uW_spline_theta(u, W, x, theta, theta_f);
n_eps = W_inte .* exp(u_inte/eps);
n_theta = n_int_x(n_eps, x, theta_f, oll);
figure(100)
h1 = plot(theta_f, n_theta, 'b-', 'LineWidth', 1.5);
hold on

figure(10)
plot(theta_f, exp(u_inte/eps),'k-.');
hold on

[max_y, idx] = max(n_theta);
max_x = theta_f(idx);
figure(100)
plot(max_x, max_y, 'bo', 'MarkerSize', 3, 'LineWidth', 2);
text(max_x, max_y, sprintf('Max: %.2f', max_y), ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','right');

theta_coarse = theta;
[~, pos] = min(abs(theta_f' - theta_coarse), [], 1);
n_theta_0 = n_theta(pos);
plot(theta_coarse, n_theta_0, '>', 'markersize', 7, 'LineWidth', 1.4);

n_nointe = W .* exp(u/eps);
n_theta_0 = (x(2) - x(1)) * (0.5 * n_nointe(1, :) + sum(n_nointe(2:end-1, :), 1) + 0.5 * n_nointe(end, :)); %#ok<NASGU>
xlabel('$\theta$', 'Interpreter', 'latex');
ylabel('$\int_0^1 n^\varepsilon dx$', 'Interpreter', 'latex');

NthDirect = 1000;
theta_o = linspace(0, 1, NthDirect + 1);
theta_o(end) = [];
marker = 1; %#ok<NASGU>
directRoot = utils.output_dir('legacy_original_bio');
directDir = utils.find_case_dir(directRoot, sprintf('eps_%s', utils.time_token(eps)), 'eps_*', true);
[directFile, directTime] = utils.find_legacy_snapshot_file(directDir, 'legacy_direct', Nx, NthDirect, t, true);
Sdir = load(directFile);
ne = Sdir.ne;
fprintf('[legacy compare] direct snapshot: %s, t=%g\n', directFile, directTime);

ne = [ne, ne(:, 1)];
netheta = n_int_x(ne, x, theta_o, oll);
figure(100)
h3 = plot(theta_o, netheta, 'r-.', 'LineWidth', 1.5);
axis([0 1 0 20])
legend([h1, h3], '$$ We^{u/\varepsilon}$$', ...
                     '$$exact ~~solution$$', 'Interpreter', 'latex');
set(gca, 'fontsize', 12);
max(netheta)
max(netheta) - max_y
hold on
title(sprintf('Nth = %d, eps = %.4g', NthDirect, eps));
