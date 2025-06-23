function evolution_compare
t = 10;   
Nth = 10;  
Nx = 10;    
ol = 1;
theta = linspace(0, 1, Nth +1); % 网格点theta
x = linspace(0, 1, Nx + 1); % 网格点x  
x = x';  
eps = 1e-2;
theta_f = 0:0.001:1; 
path = './data/';
u = load(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(ol), '.mat'));
W = load(strcat(path, 'w_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(ol), '.mat'));
u = u.u;
W = W.W; 
[u_inte, W_inte] = uW_spline_theta(u, W, x, theta, theta_f);

n_eps = W_inte .* exp(u_inte/eps);
n_theta = (x(2) - x(1)) * (0.5*n_eps(1, :) + sum(n_eps(2:end-1, :), 1) + 0.5*n_eps(end, :));
% n_theta = n_int_x(n_eps, x, theta_f);
h1 = plot(theta_f, n_theta, 'b-', 'LineWidth', 1.5);
hold on

%% 标记spline插值后的最大值点
[max_y, idx] = max(n_theta);
max_x = theta_f(idx);
% 在图上画红色圆点标记最大值位置
plot(max_x, max_y, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
% 在旁边加文字标注
text(max_x, max_y, sprintf('Max: %.2f', max_y), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');

theta_coarse = linspace(0, 1, 7 +1); % 网格点theta
[~, pos] = min(abs(theta_f'-theta_coarse), [], 1);  
n_theta_0 = n_theta(pos);
% plot(theta_coarse, n_theta_0, '>', 'markersize', 7, 'LineWidth', 1.4);
  
n_nointe = W .* exp(u/eps);   
n_theta_0 = (x(2) - x(1)) * (0.5 * n_nointe(1, :) + sum(n_nointe(2:end-1, :), 1) + 0.5 * n_nointe(end, :));
% plot(theta, n_theta_0, 'b--', 'LineWidth', 1.4);
xlabel('$\theta$', 'Interpreter', 'latex');
ylabel('$\int_0^1 n^\varepsilon dx$', 'Interpreter', 'latex');
   
path = './original_data/';
% Nth = 99;  
% Nx = 7; 
theta = linspace(0, 1, Nth +1); % 网格点theta
ne = load(strcat(path, 'ne_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'));
ne = ne.ne;   
ne = [ne, ne(:, 1)]; 
netheta = (x(2) - x(1)) * (0.5*ne(1, :) + sum(ne(2:end-1, :), 1) + 0.5*ne(end, :));
% netheta = n_int_x(ne, x, theta);
h2 = plot(theta, netheta, 'g-.', 'LineWidth', 1.5); 
axis([0 1 0 100])  
legend([h1, h2], '$$ we^{u/\varepsilon}$$', '$$ n_\varepsilon$$', 'Interpreter', 'latex');
set(gca, 'fontsize', 12);

title(sprintf('Nth =%d, eps = %.4f', Nth, eps));

end

