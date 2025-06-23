function [C, u_normal] = normalize_u(u, eps, theta)

%% u����theta��[0,1]�ϵĻ���
I = (theta(2) - theta(1)) * sum(exp(u(1:end)/eps));
%% ȷ����һ���ĳ���
C = 1 / I; 
%% ��u���й�һ��
u_normal = u + log(C) * eps;   
end

