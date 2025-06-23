function [C, u_normal] = normalize_u(u, eps, theta)

%% u关于theta在[0,1]上的积分
I = (theta(2) - theta(1)) * sum(exp(u(1:end)/eps));
%% 确定归一化的常数
C = 1 / I; 
%% 对u进行归一化
u_normal = u + log(C) * eps;   
end

