function evolution_compare
t = 50;
Nx = 10;
Nth = 10;
theta = linspace(0, 1, Nth +1); % 网格点theta
x = linspace(0, 1, Nx + 1); % 网格点x  
x = x';  
[xx, thth] = meshgrid(x, theta);
epsilon = 1e-4;
path = './data/';
load(strcat(path, 'u_', num2str(epsilon), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'));
load(strcat(path, 'w_', num2str(epsilon), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'));
u = [u, u(1)];
w = [w, w(:, 1)];
thinte = 0:0.001:1; 
xinte = 0:0.001:1;
[xxinte, ththinte] = meshgrid(xinte, thinte);

uinte = interp1(theta, u, thinte, 'spline');
winte = interp2(xx, thth, w', xxinte, ththinte, 'spline');
 
ne = winte .* exp(uinte/epsilon);
netheta = 0.001 * (0.5*ne(1, :) + sum(ne(2:end-1, :), 1) + 0.5*ne(end, :));
plot(thinte, netheta);
hold on
path = './original_data/';
load(strcat(path, 'ne_', num2str(epsilon), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'));
ne = [ne, ne(:, 1)];
netheta = 1/Nx * (0.5*ne(1, :) + sum(ne(2:end-1, :), 1) + 0.5*ne(end, :));
plot(theta, netheta);
axis([0.35 0.65 0 1500])
end

