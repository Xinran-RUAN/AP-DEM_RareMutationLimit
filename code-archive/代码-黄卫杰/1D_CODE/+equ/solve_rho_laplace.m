function rho = solve_rho_laplace(eps, theta, u, W, d_theta)

[~, k] = max(u); 

if k == 1
    km = length(u);
    kp = k + 1;
elseif k == length(u)
    km = k - 1;
    kp = 1;
else
    km = k - 1;
    kp = k + 1;
end


ddu = (u(kp)-2*u(k) + u(km))/d_theta^2;
theta_star = -d_theta/2 * (u(kp)-u(km))/(u(kp)-2*u(k)+u(km)) + theta(k);
   
lagrange1 = (theta_star - theta(k)) * (theta_star - theta(kp)) / ((theta(km) - theta(k)) *(theta(km) - theta(kp)));
lagrange2 = (theta_star - theta(km)) * (theta_star - theta(kp)) / ((theta(k) - theta(km)) *(theta(k) - theta(kp)));
lagrange3 = (theta_star - theta(km)) * (theta_star - theta(k)) / ((theta(kp) - theta(km)) *(theta(kp) - theta(k)));

W_theta_star = W(:, km) * lagrange1 + W(:, k) * lagrange2 + W(:, kp) * lagrange3;
U_theta_star = u(km) * lagrange1 + u(k) * lagrange2 + u(kp) * lagrange3;

rho = W_theta_star * exp(U_theta_star/eps) * sqrt(2*pi*eps/abs(ddu)); 
  
end
