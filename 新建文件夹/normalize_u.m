function [W_normal, u_normal] = normalize_u(W, u, eps, theta)

%% u����theta��[0,1]�ϵĻ���
% I1 = (theta(2) - theta(1)) * sum(exp(u(1:end)/eps));
%% spline ��ֵ��Ļ���
pp = spline(theta(2:end), u(2:end));
theta_f = theta(2):0.001:theta(end);
u_inte = ppval(pp, theta_f);
I = (theta_f(2)-theta_f(1)) * (0.5 * exp(u_inte(1)/eps) + ...
                               sum(exp(u_inte(2:end-1)/eps)) + ...
                               0.5 * exp(u_inte(end)/eps)) + ...
     (theta(2) - theta(1)) * 0.5 * (exp(u(1)/eps) + exp(u(2)/eps)) + ...
     (theta(2) - theta(1)) * 0.5 * (exp(u(end)/eps) + exp(u(1)/eps));
 
%% ȷ����һ���ĳ���
C = 1 / I; 
%% ��u���й�һ��
u_normal = u + log(C) * eps;  
W_normal = W / C;
end

