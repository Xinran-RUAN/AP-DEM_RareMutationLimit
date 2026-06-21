function plot_all(mesh_divide_elements, mesh_divide_nodes,...
                     theta, rho, t, u, W, H, eps, err)
    figure(1);   
    nx = 200;
    ny=200;
    [xq, yq] = meshgrid(linspace(min(mesh_divide_nodes(:, 1)), max(mesh_divide_nodes(:, 1)), nx),...
                        linspace(min(mesh_divide_nodes(:, 2)), max(mesh_divide_nodes(:, 2)), ny));
    zq = griddata(mesh_divide_nodes(:, 1), mesh_divide_nodes(:, 2), rho, xq, yq, 'natural');
    surf(xq, yq, zq);
    shading interp;  
    colorbar;
    axis tight
  
                    %  scatter3(mesh_divide_nodes(:, 1), mesh_divide_nodes(:, 2),...
  %                      rho);  
  %  plot_rho(mesh_divide_elements, mesh_divide_nodes, rho);
    xlabel('$x$', 'Interpreter', 'latex');    
    ylabel('$y$', 'Interpreter', 'latex');
    title(['Plot of $\rho(x)$, time = ', num2str(t), ', eps = ', num2str(eps)], 'Interpreter', 'latex');
    set(gcf, 'Position', [10, 50, 420, 300]);
  %  axis([0 1 0 15]); 
     dtheta_f = 1e-4;
     theta_f = 0:dtheta_f:1;
     [u_inte, W_inte] = uW_spline_theta(u, W, theta, dtheta_f);

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
%     figure(6)
%     plot([theta, 1], exp([u, u(1)]/eps));
%     set(gcf, 'Position', [10, 430, 420, 300])
    
    figure(3)
    plot(theta, H);
    hold off
    xlabel('$\theta$', 'Interpreter', 'latex');
    title(['Plot of $H(\theta)$, time = ', num2str(t), ', eps = ', num2str(eps)], 'Interpreter', 'latex');
    set(gcf, 'Position', [870, 50, 420, 300])
   % axis([0 1 -1 2]);
    
    figure(5)
    node_num = length(mesh_divide_nodes);
    plot([theta, 1], [W(ceil(node_num/2), :), W(ceil(node_num/2), 1)]);
    hold off
    xlabel('$\theta$', 'Interpreter', 'latex');
    title(['Plot of $W(\theta)$, time = ', num2str(t), ', eps = ', num2str(eps)], 'Interpreter', 'latex');
    set(gcf, 'Position', [440, 430, 420, 300])
      
    %% »¹Ô­ne
    figure(4);    
    hold off
    ne = W .* exp(u/eps);   
    Nth = length(u); 
    ntheta = int_theta_n(mesh_divide_elements, mesh_divide_nodes, ne, Nth);
    plot([theta,1], [ntheta, ntheta(1)]); 
    xlabel('$\theta$', 'Interpreter', 'latex');
    title(['Plot of $\int_0^1 n^\varepsilon dx$, time = ', num2str(t), ', eps = ', num2str(eps)], 'Interpreter', 'latex');
    set(gcf, 'Position', [870, 430, 420, 300])
    hold on

     n_f = W_inte .* exp(u_inte/eps);   
     Nth_f = length(u_inte);
     ntheta_f = int_theta_n(mesh_divide_elements, mesh_divide_nodes, n_f, Nth_f);
     plot(theta_f, ntheta_f, 'r-'); 
end

