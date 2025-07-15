function u_new = solve_u(u, u_temp, B, H, d_theta, dt, theta, ol)
%% �����ʽȫ������������u��ʾ��һ��ʱ�䲽��ֵ��u_temp��ʾ���������е�ֵ
% �������������������Ϣp^2����������һ�׵�p>0,������/�����֣������Ҳ��
grad_sq = weno5_h_hat(u, d_theta);
  
% ���׵�����ʽ��һ�׵���ӭ��,��ʽEuler
u_new = B\(u + dt * (- grad_sq - H))';

u_new = u_new';
end

