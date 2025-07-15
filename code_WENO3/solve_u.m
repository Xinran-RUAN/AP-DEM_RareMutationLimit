function u_new = solve_u(u, u_temp, B, H, d_theta, dt, theta, ol)
%% 如果格式全隐，需迭代求解u表示上一个时间步的值，u_temp表示迭代过程中的值
% 非线性项的正负处理：信息p^2从左来，即一阶导p>0,则左差分/后向差分；否则右差分
grad_sq = weno5_h_hat(u, d_theta);
  
% 二阶导用隐式，一阶导用迎风,隐式Euler
u_new = B\(u + dt * (- grad_sq - H))';

u_new = u_new';
end

