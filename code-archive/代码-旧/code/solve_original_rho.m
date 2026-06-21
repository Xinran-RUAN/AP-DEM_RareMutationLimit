<<<<<<< HEAD
function rho = solve_original_rho(ne, Nth, x, theta)
%% 插值准备
ne_aux = [ne, ne(:, 1)];
%% spline 插值
[xx, thth] = meshgrid(x, [theta, 1]);
th_f = 0:0.0001:1; 
[X_f, Th_f] = meshgrid(x, th_f);

ne_f = interp2(xx, thth, ne_aux', X_f, Th_f, 'spline');
ne_f = ne_f'; 

rho = (Th_f(2) - Th_f(1)) * sum(ne_f(:, 1: end-1), 2);

% rho = (theta(2) - theta(1)) * sum(ne(:, 1: end-1), 2); 
end 

=======
function rho = solve_original_rho(ne, Nth, x, theta)
%% 插值准备
ne_aux = [ne, ne(:, 1)];
%% spline 插值
[xx, thth] = meshgrid(x, [theta, 1]);
th_f = 0:0.0001:1; 
[X_f, Th_f] = meshgrid(x, th_f);

ne_f = interp2(xx, thth, ne_aux', X_f, Th_f, 'spline');
ne_f = ne_f'; 

rho = (Th_f(2) - Th_f(1)) * sum(ne_f(:, 1: end-1), 2);

% rho = (theta(2) - theta(1)) * sum(ne(:, 1: end-1), 2); 
end 

>>>>>>> d725b20c4e0dc455060136e6e0d3a79bba8f525a
