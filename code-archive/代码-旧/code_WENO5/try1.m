function try1
N = 40;
x = linspace(0,1,N);
dx = x(2)-x(1);
x(end) = [];
 u = sin(2*pi*x);  
du_exact = 2*pi*cos(2*pi*x);

% u = abs(x-0.5);

pL = WENO5_left_RXR(u,dx);
pR = WENO5_right_RXR(u, dx);

plot(x, pL - du_exact);
 
figure(10)
du_num = RF_flux(pR, pL);
 
% du_num = weno3_diff(u, dx);

du_num2 = RF_flux_neg_uv(pR, pL, pR, pL);

figure(100)
plot(x, du_num2+du_exact.^2);
hold on  
plot(x, du_num, 'r--');
plot(x, du_num2, 'k-');
% plot(x, du_num2, 'k-.');
% max(abs(du_num -du_num2))
% figure(200)
% plot(x, du_num -du_num2);


error = max(abs(du_num - du_exact));
fprintf('×î´óÎó²î: %.3e\n', error)

end

