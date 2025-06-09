function[grad_sq] = compute_grad_sq(u_temp, d_theta)
u_plus = [u_temp(2:end), u_temp(1)]; % 关于theta是周期边界
u_minus = [u_temp(end), u_temp(1:end-1)];
% 非线性项的正负处理：信息p^2从左来，即一阶导p>0,则左差分/后向差分；否则右差分
pL = (u_temp - u_minus) / d_theta; % 左差分 u_j - u_{j-1}
pR = (u_plus - u_temp) / d_theta; % 右差分 u_{j+1} - u_j
grad_sq = min(pL.^2, pR.^2) .* (pL .* pR >= 0) + ...
          pL .* 0 .* (pL .* pR < 0);
grad_sq(1) = pR(1).^2;
