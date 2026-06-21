function plot_all(x, theta, rho, t, dx, u, W, H, eps, N_x, err)
figure(1);              
    plot(x, rho);  
    xlabel('$x$', 'Interpreter', 'latex');    
    title(['Plot of $\rho(x)$, time = ', num2str(t), ', eps = ', num2str(eps)], 'Interpreter', 'latex');
    set(gcf, 'Position', [10, 50, 420, 300]);
  %  axis([0 1 0 15]);
  eps = 1e-2;
  theta_f = 0:0.001:1;
  [u_inte, W_inte] = uW_spline_theta(u, W, x, [theta, 1], theta_f);

    %% »­Í¼  
    figure(2); 
    plot([theta, 1], [u, u(1)]); 
    hold on
    plot(theta_f, u_inte);
    hold off
    xlabel('$\theta$', 'Interpreter', 'latex');
    title(['Plot of $u(\theta)$, time = ', num2str(t), ', eps = ', num2str(eps)], 'Interpreter', 'latex');
    set(gcf, 'Position', [440, 50, 420, 300])
    % axis([0 1 -0.5 0.5]);
    
    figure(3)
    plot(theta, H);
    xlabel('$\theta$', 'Interpreter', 'latex');
    title(['Plot of $H(\theta)$, time = ', num2str(t), ', eps = ', num2str(eps)], 'Interpreter', 'latex');
    set(gcf, 'Position', [870, 50, 420, 300])
   % axis([0 1 -1 2]);
    
    figure(5)
    plot([theta, 1], [W(ceil(N_x/2), :), W(ceil(N_x/2), 1)]);
    hold on
    plot(theta_f, W_inte(ceil(N_x/2), :));
    hold off
    xlabel('$\theta$', 'Interpreter', 'latex');
    title(['Plot of $W(\theta)$, time = ', num2str(t), ', eps = ', num2str(eps)], 'Interpreter', 'latex');
    set(gcf, 'Position', [440, 430, 420, 300])
      
    %% »¹Ô­ne
    figure(4);    
    hold off
    ne = W .* exp(u/eps);    
    ntheta = dx * (ne(1, :)/2 + sum(ne(2:N_x, :), 1) + ne(end, :)/2);
    plot([theta,1], [ntheta, ntheta(1)]); 
    xlabel('$\theta$', 'Interpreter', 'latex');
    title(['Plot of $\int_0^1 n^\varepsilon dx$, time = ', num2str(t), ', eps = ', num2str(eps)], 'Interpreter', 'latex');
    set(gcf, 'Position', [870, 430, 420, 300])
    hold on
    d_theta_f = 1e-3;
    dx_f = 1e-3;
    [X, Theta] = meshgrid(x, [theta, 1]);
    theta_f = 0:d_theta_f:1;
    x_f = 0: dx_f:1;
    [X_f, Theta_f] = meshgrid(x_f, theta_f);
    u_aux = [u, u(:, 1)];
    w_aux = [W, W(:, 1)];  
    uinte = interp1([theta, 1], u_aux, theta_f, 'spline');
    winte = interp2(X, Theta, w_aux', X_f, Theta_f, 'spline');
    winte = winte';
    ne = winte .* exp(uinte/eps);    
    ntheta = dx_f * (ne(1, :)/2 + sum(ne(1:end-1, :), 1) + ne(end, :)/2);
    plot(Theta_f, ntheta, 'r-'); 
end

