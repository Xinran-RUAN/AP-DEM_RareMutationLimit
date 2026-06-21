function [u_f, W_f] = uW_spline_theta(u, W, theta, dtheta_f)

     theta_f = 0:dtheta_f:1;
     theta_ext = [theta, 1];
     u_ext = [u, u(1)];
     u_f = interp1(theta_ext, u_ext, theta_f, 'spline');
     
     %% WµÄspline²åÖµ
     W_ext = [W, W(:, 1)];
     Wf_tmp = interp1(theta_ext, W_ext', theta_f, 'spline');
     W_f = Wf_tmp';  

end

