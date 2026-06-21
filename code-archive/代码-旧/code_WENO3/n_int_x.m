function n_theta = n_int_x(n, x, theta_f)
%% spline ²åÖµ
[xx, thth] = meshgrid(x, theta_f);
xinte = 0:0.001:1; 
[xxinte, ththinte] = meshgrid(xinte, theta_f);
   
neinte = interp2(xx, thth, n', xxinte, ththinte, 'spline');
neinte = neinte'; 
   
n_theta = (xinte(2) - xinte(1)) * (neinte(1, :) / 2 + ...
                                   sum(neinte(2:end-1, :), 1) + ...
                                   neinte(end, :) / 2);
end

