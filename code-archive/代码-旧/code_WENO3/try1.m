function try1
N = 100;
x = linspace(0,1,N);
dx = x(2)-x(1);
x(end) = [];
u = sin(2*pi*x);  
du_exact = 2*pi*cos(2*pi*x);

pL = weno3_left(u,dx);
pR = weno3_right(u, dx);
du_num = RF_flux(pR, pL);

du_num2 = RF_flux_neg_uv(pR, pL, pR, pL);

figure(100)
plot(x, du_exact);
hold on
plot(x, pL, 'r--');
plot(x, pR, 'k-.');
% plot(x, du_num2, 'k-.');
% max(abs(du_num -du_num2))
% figure(200)
% plot(x, du_num -du_num2);


error = max(abs(du_num - du_exact));
fprintf('×î´óÎó²î: %.3e\n', error)

end

