function rho = solve_original_rho(ne, x, theta)
%% 插值准备
ne_aux = [ne, ne(:, 1)];
%% spline 插值
[xx, thth] = meshgrid(x, [theta, 1]);
th_f = 0:0.0001:1; 
[X_f, Th_f] = meshgrid(x, th_f);

ne_f = interp2(xx, thth, ne_aux', X_f, Th_f, 'spline');
ne_f = ne_f'; 

rho = (Th_f(2) - Th_f(1)) * sum(ne_f(:, 1: end-1), 2);
end 

