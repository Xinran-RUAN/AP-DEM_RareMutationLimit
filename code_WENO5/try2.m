function try2
% ����
Nx = 1000; 
x = linspace(0,1,Nx);
dx = x(2)-x(1);


% ����
u = sin(2*pi*x);
v = cos(4 *pi*x);

% �浼��
u_x_exact = 2*pi*cos(2*pi*x);
v_x_exact = - 4 *pi*sin(4 * pi*x);

% ��ֵ�������򵥲�֣���ɻ��� WENO ���ർ����
u_x_minus = weno5_left(u, dx);
u_x_plus = weno5_right(u, dx);
v_x_minus = weno5_left(v, dx);
v_x_plus = weno5_right(v, dx);

% ��ֵ flux
Hhat = RF_flux_neg_uv(u_x_plus, u_x_minus, v_x_plus, v_x_minus);

% �� flux
Htrue = -u_x_exact .* v_x_exact;

plot(x, Htrue);
hold on
plot(x, Hhat, 'r--');
end