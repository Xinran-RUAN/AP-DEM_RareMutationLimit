function evolution_compare
close all
t = 5;
Nx = 7;  
Nth = 7;
theta = linspace(0, 1, Nth +1); % 网格点theta
x = linspace(0, 1, Nx + 1); % 网格点x  
x = x';  
[X, Theta] = meshgrid(x, theta);
eps = 1e-2;
path = './data/';
u = load(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'));
W = load(strcat(path, 'w_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'));
u = u.u;
W = W.W;
u = [u, u(1)];
W = [W, W(:, 1)];
theta_f = 0:0.001:1; 
[X_f, Theta_f] = meshgrid(x, theta_f);

u_inte = interp1(theta, u, theta_f, 'spline');
W_inte = interp2(X, Theta, W', X_f, Theta_f, 'spline');
 
n_eps = W_inte' .* exp(u_inte/eps);
n_theta = (x(2) - x(1)) * (0.5*n_eps(1, :) + sum(n_eps(2:end-1, :), 1) + 0.5*n_eps(end, :));
plot(theta_f, n_theta);
hold on

[~, pos] = min(abs(theta_f'-theta), [], 1);  
n_theta_0 = n_theta(pos);
plot(theta, n_theta_0, '>');
  
n_nointe = W .* exp(u/eps);   
n_theta_0 = (x(2) - x(1)) * (0.5 * n_nointe(1, :) + sum(n_nointe(2:end-1, :), 1) + 0.5 * n_nointe(end, :));
plot(theta, n_theta_0, '--');
xlabel('$\theta$', 'Interpreter', 'latex');
ylabel('$\int_0^1 n^\varepsilon dx$', 'Interpreter', 'latex');

path = './original_data/';
t=49; 
ne = load(strcat(path, 'ne_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'));
ne = ne.ne; 
ne = [ne, ne(:, 1)];     
netheta = 1/Nx * (0.5*ne(1, :) + sum(ne(2:end-1, :), 1) + 0.5*ne(end, :));
plot(theta, netheta); 
axis([0 1 0 100])    
end

