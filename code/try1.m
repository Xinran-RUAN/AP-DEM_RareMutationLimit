function try1
N = 100;
x = linspace(0,1,N);
dx = x(2)-x(1);
x(end) = [];
u = sin(2*pi*x);  
du_exact = 2*pi*cos(2*pi*x);

pL = weno5_left(u,dx);
pR = weno5_right(u, dx);
du_num = RF_flux(pR, pL);

plot(x, -du_exact.^2);
hold on
plot(x, du_num, 'r--');
error = max(abs(du_num - du_exact));
fprintf('×î´óÎó²î: %.3e\n', error)

end

