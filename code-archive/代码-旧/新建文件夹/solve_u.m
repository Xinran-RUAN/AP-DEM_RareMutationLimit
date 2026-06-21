function u_new = solve_u(u, u_temp, B, H, d_theta, dt, theta, ol)
%% 如果格式全隐，需迭代求解u表示上一个时间步的值，u_temp表示迭代过程中的值
u_plus = [u_temp(2:end), u_temp(1)]; % 关于theta是周期边界
u_minus = [u_temp(end), u_temp(1:end-1)];
% 非线性项的正负处理：信息p^2从左来，即一阶导p>0,则左差分/后向差分；否则右差分

if ol == 0
    %% 取绝对值小的导数
    pL = (u_temp - u_minus) / d_theta; % 左差分 u_j - u_{j-1}
    pR = (u_plus - u_temp) / d_theta; % 右差分 u_{j+1} - u_j
    grad_sq = min(pL.^2, pR.^2) .* (pL .* pR >= 0) + ...
    pL .* 0 .* (pL .* pR < 0);
    grad_sq(1) = pR(1).^2;
%     %% WENO5
      du = weno5_diff(u, d_theta);
      grad_sq = du.^2;
else
    if ol == 1
        %% spline插值求导数
        du = 0 .* u;
        pp = spline(theta(2:end), u(2:end));
        pp_der = fnder(pp, 1);
        du(2:end) = ppval(pp_der, theta(2:end));
        pR = (u_plus - u_temp) / d_theta; % 右差分 u_{j+1} - u_j
       %  du(1) = pR(1);
        grad_sq = du.^2;
    else
        pp = spline([theta, 1], [u, u(1)]);
        pp_der = fnder(pp, 1);
        du = ppval(pp_der, theta);
        grad_sq = du.^2;
    end
end
  
% 二阶导用隐式，一阶导用迎风,隐式Euler
u_new = B\(u + dt * (grad_sq - H))';

u_new = u_new';
end

