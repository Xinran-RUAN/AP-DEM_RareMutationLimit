<<<<<<< HEAD
function u_new = solve_u_semi(u, u_temp, B, H, d_theta, dt, theta, ol)
%% �����ʽȫ������������u��ʾ��һ��ʱ�䲽��ֵ��u_temp��ʾ���������е�ֵ
u_plus = [u_temp(2:end), u_temp(1)]; % ����theta�����ڱ߽�
u_minus = [u_temp(end), u_temp(1:end-1)];
N_theta = length(u);  
% �������������������Ϣp^2����������һ�׵�p>0,������/�����֣������Ҳ��
if ol == 0
    %% ȡ����ֵС�ĵ���
    pL = (u_temp - u_minus) / d_theta; % ���� u_j - u_{j-1}
    pR = (u_plus - u_temp) / d_theta; % �Ҳ�� u_{j+1} - u_j
    du = min(pL, pR) .* (pL .* pR > 0) .* (pL > 0) + ...
         max(pL, pR) .* (pL .* pR > 0) .* (pL < 0) + ...
         pL .* 0 .* (pL .* pR <= 0);
    du(1) = pR(1);
else
    if ol == 1
        %% spline��ֵ����
        du = 0 .* u;
        pp = spline(theta(2:end), u(2:end));
        pp_der = fnder(pp, 1);
        du(2:end) = ppval(pp_der, theta(2:end));
        pR = (u_plus - u_temp) / d_theta; % �Ҳ�� u_{j+1} - u_j
        du(1) = pR(1);
    else
        pp = spline([theta, 1], [u, u(1)]);
        pp_der = fnder(pp, 1);
        du = ppval(pp_der, theta);
    end
end
% ���׵�����ʽ��һ�׵���ӭ��,��ʽEuler 
B_semi_1d = du' .* (dt / (2.*d_theta)) .* (- diag(ones(N_theta-1, 1), 1) +...
                                  + diag(ones(N_theta-1, 1), -1));
B_semi_1d(1, N_theta) = dt / (2 * d_theta) .* du(1);
B_semi_1d(N_theta, 1) = - dt / (2 * d_theta) .* du(end);  
% B_semi_1d(1, 1) = dt / (2 * d_theta) .* du(1);
% B_semi_1d(N_theta, N_theta) = - dt / (2 * d_theta) .* du(end);  

u_new = (B + B_semi_1d)\(u + dt * (- H))';

u_new = u_new';
end  

=======
function u_new = solve_u_semi(u, u_temp, B, H, d_theta, dt, theta, ol)
%% �����ʽȫ������������u��ʾ��һ��ʱ�䲽��ֵ��u_temp��ʾ���������е�ֵ
u_plus = [u_temp(2:end), u_temp(1)]; % ����theta�����ڱ߽�
u_minus = [u_temp(end), u_temp(1:end-1)];
N_theta = length(u);  
% �������������������Ϣp^2����������һ�׵�p>0,������/�����֣������Ҳ��
if ol == 0
    %% ȡ����ֵС�ĵ���
    pL = (u_temp - u_minus) / d_theta; % ���� u_j - u_{j-1}
    pR = (u_plus - u_temp) / d_theta; % �Ҳ�� u_{j+1} - u_j
    du = min(pL, pR) .* (pL .* pR > 0) .* (pL > 0) + ...
         max(pL, pR) .* (pL .* pR > 0) .* (pL < 0) + ...
         pL .* 0 .* (pL .* pR <= 0);
    du(1) = pR(1);
else
    if ol == 1
        %% spline��ֵ����
        du = 0 .* u;
        pp = spline(theta(2:end), u(2:end));
        pp_der = fnder(pp, 1);
        du(2:end) = ppval(pp_der, theta(2:end));
        pR = (u_plus - u_temp) / d_theta; % �Ҳ�� u_{j+1} - u_j
        du(1) = pR(1);
    else
        pp = spline([theta, 1], [u, u(1)]);
        pp_der = fnder(pp, 1);
        du = ppval(pp_der, theta);
    end
end
% ���׵�����ʽ��һ�׵���ӭ��,��ʽEuler 
B_semi_1d = du' .* (dt / (2.*d_theta)) .* (- diag(ones(N_theta-1, 1), 1) +...
                                  + diag(ones(N_theta-1, 1), -1));
B_semi_1d(1, N_theta) = dt / (2 * d_theta) .* du(1);
B_semi_1d(N_theta, 1) = - dt / (2 * d_theta) .* du(end);  
% B_semi_1d(1, 1) = dt / (2 * d_theta) .* du(1);
% B_semi_1d(N_theta, N_theta) = - dt / (2 * d_theta) .* du(end);  

u_new = (B + B_semi_1d)\(u + dt * (- H))';

u_new = u_new';
end  

>>>>>>> d725b20c4e0dc455060136e6e0d3a79bba8f525a
