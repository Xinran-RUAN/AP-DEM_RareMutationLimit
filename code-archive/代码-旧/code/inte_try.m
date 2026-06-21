function inte_try

dth = [0.02, 0.01, 0.001, 0.0001, 0.00001];
rho = zeros(1, length(dth));
for kk = 1: length(dth)
theta = 0:dth(kk):1;
epsilon = 1e-5;
n = exp((sin(pi*theta) -1) / epsilon);
rho(kk) = (theta(2)-theta(1)) * (0.5*n(1) + sum(n(2:end-1), 2) + 0.5*n(end));
end

% plot(dth, rho, '-*');

figure(10)
plot(theta, n, 'b-');
hold on
axis([0.48 0.52 0 1]);

dtheta = 0.01;
theta = 0:dtheta:1;
u = sin(pi*theta) -1;
n = exp(u / epsilon);
figure(10)
plot(theta, n, 'r--');
hold on
thth = 0:0.00001:1;
uu = interp1(theta, u, thth, 'spline');
nn = exp(uu / epsilon); 
figure(10)
plot(thth, nn, 'g-.');
hold on

rhoo = (thth(2)-thth(1)) * (0.5*nn(1) + sum(nn(2:end-1), 2) + 0.5*nn(end));
hold on 
plot(thth(2)-thth(1), rhoo, 'r>'); 
end

