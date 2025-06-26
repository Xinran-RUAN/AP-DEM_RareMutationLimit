function u_new = solve_u(u, u_temp, B, H, d_theta, dt, is_symmetric_B)
grad_sq = compute_grad_sq(u_temp, d_theta);

%    grad_sq1 = (Du(u) / d_theta).^2;
%     pp = spline([theta, 1], [u, u(1)]);
%     pp_der = fnder(pp, 1);
%     du = ppval(pp_der, theta);
%     grad_sq = du.^2;
% 二阶导用隐式，一阶导用迎风,隐式Euler
RHS = u + dt * (grad_sq - H);
% if is_symmetric_B
%     try 
%         R = chol(B, 'lower');
%         u_new = R' \ (R \ RHS');
%     catch
%         warning('Cholesky failed, fallback to LU');
%         u_new = B \ RHS';
%     end
% else
%     u_new = B \ RHS';
% end
u_new = B \ RHS';
%% fixed iteration for grad_sq
% tol = 1e-14;
% err = 1;
% u_o = u_temp;
% while err > tol
%     grad_sq_o = compute_grad_sq(u_o, d_theta);
%     RHS_o = u_temp + dt * (grad_sq_o - H);
%     u_new_predict = B \ RHS_o';
%     u_t = u_new_predict';
%     grad_sq_n = compute_grad_sq(u_t, d_theta);
%     RHS_n = u_temp + dt * (grad_sq_n - H);
%     u_new = B \ RHS_n';
%     u_o = u_new';
%     err = max(abs(u_o - u_t));
% end

%% for test - artificial u
% N_theta = 7; % theta方向网格步长1/Nth
% d_theta = 1 / N_theta;      
% theta = 0:d_theta:1-d_theta;
% H = (theta - 0.5).^2;
% H(2:4) = H(end:-1:end-2);
% grad_sq(2:4) = grad_sq(end:-1:end-2);
% u_new = B\(u + dt * (grad_sq - H))';

%%
% plot(u_new(2:4) - u_new(end:-1:end-2)); hold on
% plot(H(2:4) - H(end:-1:end-2),'r--'); hold off
% pause(0.001)

% if abs(u_new(end-1)-u_new(3)) > 1e-6
%     pause(0.1)
% end

u_new = u_new';
end

