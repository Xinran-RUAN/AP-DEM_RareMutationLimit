function u_new = solve_u(u, u_temp, B, H, d_theta, dt, theta, ol)
%% 如果格式全隐，需迭代求解u表示上一个时间步的值，u_temp表示迭代过程中的值
u_plus = [u_temp(2:end), u_temp(1)]; 
u_minus = [u_temp(end), u_temp(1:end-1)];
% 非线性项的正负处理：信息p^2从左来，即一阶导p>0,则左差分/后向差分；否则右差分

%% spline插值求导数
pp = spline(theta(2:end-1), u(2:end-1));
pp_der = fnder(pp, 1);
du = ppval(pp_der, theta(2:end-1));
grad_sq = du.^2;
   
% 二阶导用隐式，一阶导用迎风,隐式Euler
u_new = B\(u(2:end-1) + dt * (grad_sq - H(2:end)))';

u_new = [u_new(1); u_new; u_new(end)]';
end

