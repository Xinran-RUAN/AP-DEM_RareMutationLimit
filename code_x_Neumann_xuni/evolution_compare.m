function evolution_compare
t = 2;  
eps = 1e-2; 
ol = 0; 
oll = 1;

%% �������
Nth = 15;      
Nx = 10;         
theta = linspace(0, 1, Nth +1); % �����theta
x = linspace(0, 1, Nx + 1); % �����x  
x = x';  
theta_f = 0:0.001:1; 
%% W,u ����
path = './data/';
u = load(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(ol), '.mat'));
W = load(strcat(path, 'w_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(ol), '.mat'));
u = u.u;  
W = W.W;  

%% spline��ֵ  
[u_inte, W_inte] = uW_spline_theta(u, W, x, theta, theta_f);
n_eps = W_inte .* exp(u_inte/eps);
n_theta = n_int_x(n_eps, x, theta_f, oll);
h1 = plot(theta_f, n_theta, 'b-', 'LineWidth', 1.5);
hold on

%% ���spline��ֵ������ֵ��
[max_y, idx] = max(n_theta);
max_x = theta_f(idx);
% ��ͼ�ϻ���ɫԲ�������ֵλ��
plot(max_x, max_y, 'bo', 'MarkerSize', 3, 'LineWidth', 2);
% ���Ա߼����ֱ�ע
text(max_x, max_y, sprintf('Max: %.2f', max_y), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');

%% ��ɢ�㻭һ��
theta_coarse = theta;%linspace(0, 1, 7 +1); % �����theta
[~, pos] = min(abs(theta_f'-theta_coarse), [], 1);  
n_theta_0 = n_theta(pos);
plot(theta_coarse, n_theta_0, '>', 'markersize', 7, 'LineWidth', 1.4);
  
n_nointe = W .* exp(u/eps);   
n_theta_0 = (x(2) - x(1)) * (0.5 * n_nointe(1, :) + sum(n_nointe(2:end-1, :), 1) + 0.5 * n_nointe(end, :));
plot(theta_coarse, [n_theta_0, n_theta_0(1)], 'b--', 'LineWidth', 1.4);
xlabel('$\theta$', 'Interpreter', 'latex');
ylabel('$\int_0^1 n^\varepsilon dx$', 'Interpreter', 'latex');

%% ԭʼ����Ľ�
path = './original_data_nomid/';
t = 5;
Nth = 1000;      
Nx = 10;  
theta_o = linspace(0, 1, Nth +1); % �����theta
x = linspace(0, 1, Nx + 1); % �����x  
x = x'; 
marker = 1;
ne = load(strcat(path, 'ne_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(marker), '.mat'));
ne = ne.ne; 
ne = [ne, ne(:, 1)];  
netheta = n_int_x(ne, x, theta_o, oll);
h3 = plot(theta_o, netheta, 'r-.', 'LineWidth', 1.5); 
axis([0 1 0 100])  
legend([h1, h3], '$$ We^{u/\varepsilon}$$',...
                     '$$exact ~~solution$$', 'Interpreter', 'latex');
set(gca, 'fontsize', 12);
max(netheta)
max(netheta) - max_y
hold on
% max_y

title(sprintf('Nth =%d, eps = %.4f', Nth, eps));

end