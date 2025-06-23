function rho = solve_original_rho(ne, Nth, x, theta)
%% ��ֵ׼��
ne_aux = [ne, ne(:, 1)];
%% spline ��ֵ
[xx, thth] = meshgrid(x, [theta, 1]);
thinte = 0:0.001:1; 
[xxinte, ththinte] = meshgrid(x, thinte);

neinte = interp2(xx, thth, ne_aux', xxinte, ththinte, 'spline');
neinte = neinte'; 
   
rho = (ththinte(2) - ththinte(1)) * sum(neinte(:, 1: end-1), 2);
end

