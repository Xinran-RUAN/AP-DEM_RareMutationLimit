function evolution_compare
t = 50;   
Nth = 10;    
Nx = 10;       
ol = 0;  
marker = 0;
theta = linspace(0, 1, Nth +1); % 网格点theta
x = linspace(0, 1, Nx + 1); % 网格点x  
x = x';  
eps = 1e-2;
theta_f = 0:0.001:1; 
path = './data/';
if marker == -1
    u = load(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(ol), '.mat'));
    W = load(strcat(path, 'w_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(ol), '.mat'));
else
    u = load(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(ol), '_', num2str(marker), '.mat'));
    W = load(strcat(path, 'w_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(ol), '_', num2str(marker), '.mat'));
end
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
plot(max_x, max_y, 'bo', 'MarkerSize', 3, 'LineWidth', 2);
% 在旁边加文字标注
text(max_x, max_y, sprintf('Max: %.2f', max_y), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');

theta_coarse = linspace(0, 1, 7 +1); % 网格点theta
[~, pos] = min(abs(theta_f'-theta_coarse), [], 1);  
n_theta_0 = n_theta(pos);
plot(theta_coarse, n_theta_0, '>', 'markersize', 7, 'LineWidth', 1.4);
  
n_nointe = [W, W(:, 1)] .* exp(u/eps);   
n_theta_0 = (x(2) - x(1)) * (0.5 * n_nointe(1, :) + sum(n_nointe(2:end-1, :), 1) + 0.5 * n_nointe(end, :));
%plot(theta, n_theta_0, 'b--', 'LineWidth', 1.4);
xlabel('$\theta$', 'Interpreter', 'latex');
ylabel('$\int_0^1 n^\varepsilon dx$', 'Interpreter', 'latex');

%% 原始问题的解
path = './original_data_nomid/';
t = 3;
Nth = 1000;      
Nx = 10; 
theta = linspace(0, 1, Nth +1); % 网格点theta
x = linspace(0, 1, Nx + 1); % 网格点x  
x = x'; 
% if marker == 0
%     ne = load(strcat(path, 'ne_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'));
% else
    marker = 1;
    ne = load(strcat(path, 'ne_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(marker), '.mat'));
% end
ne = ne.ne;
ne = [ne, ne(:, 1)]; 
netheta = (x(2) - x(1)) * (0.5*ne(1, :) + sum(ne(2:end-1, :), 1) + 0.5*ne(end, :));
h3 = plot(theta, netheta, 'g-.', 'LineWidth', 1.5); 
axis([0 1 0 100])  
legend([h1, h3], '$$ we^{u/\varepsilon}, spline~ (0,1)$$',...
                     '$$ n_\varepsilon$$', 'Interpreter', 'latex');
set(gca, 'fontsize', 12);
max(netheta)
hold on

title(sprintf('Nth =%d, eps = %.4f', Nth, eps));
%axis([0.47 0.53 60 65])  
end

