function [W_normal, u_normal] = normalize_u(W, u, eps, theta)
%% spline 插值后的积分
pp = spline(theta(2:end-1), u(2:end-1));
theta_f = theta(2):0.001:theta(end-1);
u_inte = ppval(pp, theta_f);
I = (theta_f(2)-theta_f(1)) * (0.5 * exp(u_inte(1)/eps) + ...
                               sum(exp(u_inte(2:end-1)/eps)) + ...
                               0.5 * exp(u_inte(end)/eps)) + ...
     (theta(2) - theta(1)) * 0.5 * (exp(u(1)/eps) + exp(u(2)/eps)) + ...
     (theta(2) - theta(1)) * 0.5 * (exp(u(end)/eps) + exp(u(end-1)/eps));
 
%% 确定归一化的常数
C = 1 / I; 
%% 对u进行归一化
u_normal = u + log(C) * eps;  
W_normal = W / C;
end

