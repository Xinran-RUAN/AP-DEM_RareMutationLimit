close all
clear 

t = 2;  
eps = 1e-0; 
ol = 0; 
oll = 0;

%% 网格参数
Nth = [15, 20, 50, 100, 200, 500, 1000];     
Nx = 5; 
x = linspace(0, 1, Nx + 1); % 网格点x  
x = x';  
M = zeros(1, length(Nth));
for kk = 1:length(Nth)
    
    theta = linspace(0, 1, Nth(kk) +1); % 网格点theta
    theta_f = 0:0.0001:1;
    %% W,u 数据
    path = './data/';
    u = load(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth(kk)), '_', num2str(ol), '.mat'));
    W = load(strcat(path, 'w_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth(kk)), '_', num2str(ol), '.mat'));
    u = u.u;
    W = W.W;

    %% spline插值
    [u_inte, W_inte] = uW_spline_theta(u, W, x, theta, theta_f);
%     u_inte = u;
%     W_inte = W;
    n_eps = W_inte .* exp(u_inte/eps);
    n_theta = n_int_x(n_eps, x, theta_f, oll);
   
    M(kk) = max(n_theta);
end

No = [20, 50, 100, 200, 500, 1000];
path = './original_data/';
Mo = zeros(1, length(No));

for k = 1: length(No)
    N = No(k);
    theta_o = linspace(0, 1, N +1); % 网格点theta
    marker = 1;
    ne = load(strcat(path, 'ne_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(N), '_', num2str(marker), '.mat'));
    ne = ne.ne;
    ne = [ne, ne(:, 1)];
    netheta = n_int_x(ne, x, theta_o, oll);
    Mo(k) = max(netheta);
end
 
figure
plot(log10(Nth(1:end-1)), log10(M(1:end-1)-M(end)), 'b-*');
hold on
plot(log10(No), log10(Mo-M(end)), 'r-->');

figure
plot(log10(Nth), log10(M-Mo(end)), 'b-*');
hold on
plot(log10(No(1:end-1)), log10(Mo(1:end-1)-Mo(end)), 'r-->');

