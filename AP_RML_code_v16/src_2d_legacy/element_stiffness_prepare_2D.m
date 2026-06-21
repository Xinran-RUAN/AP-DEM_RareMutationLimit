function [k1, k2_x, k2_y, con] = element_stiffness_prepare_2D(xi, eta, mesh_divide_elements, mesh_divide_nodes, N)
k1 = zeros(3, 3, length(xi));
k2_x = zeros(3, 3, length(xi), N);
k2_y = zeros(3, 3, length(xi), N);

con = zeros(3, 1, length(xi));

for e = 1: N
    nodes = mesh_divide_elements(e, :);
    x_i = mesh_divide_nodes(nodes(1), 1);
    y_i = mesh_divide_nodes(nodes(1), 2);
    x_j = mesh_divide_nodes(nodes(2), 1);
    y_j = mesh_divide_nodes(nodes(2), 2);
    x_k = mesh_divide_nodes(nodes(3), 1);
    y_k = mesh_divide_nodes(nodes(3), 2);
    
    Jacobi = (x_j - x_i) * (y_k - y_i) - (x_k - x_i) * (y_j - y_i);
    
        xi_x = (y_k - y_i) / Jacobi;
    xi_y = (x_i - x_k) / Jacobi;
    
    eta_x = (y_i - y_j) / Jacobi;
    eta_y = (x_j - x_i) / Jacobi;
    
    for l_1 = 1:3
        for l_2 = 1:3
            k2_x(l_2, l_1, :, e) = (d_varphi_1_2D(xi, eta, l_1, 1) .* xi_x + ...
                d_varphi_1_2D(xi, eta, l_1, 2) .* eta_x) .* (d_varphi_1_2D(xi, eta, l_2, 1) .* xi_x + ...
                + d_varphi_1_2D(xi, eta, l_2, 2) .* eta_x);
            k2_y(l_2, l_1, :, e) = (d_varphi_1_2D(xi, eta, l_1, 1) .* xi_y + ...
                d_varphi_1_2D(xi, eta, l_1, 2) .* eta_y) .* (d_varphi_1_2D(xi, eta, l_2, 1) .* xi_y + ...
                + d_varphi_1_2D(xi, eta, l_2, 2) .* eta_y);            
        end
    end
end

for l_1 = 1:3
    for l_2 = 1:3
        k1(l_2, l_1, :) = varphi_1_2D(xi, eta, l_1) .* varphi_1_2D(xi, eta, l_2);
    end
end

for l_1 = 1:3
    con(l_1, 1, :) = varphi_1_2D(xi, eta, l_1);
end

end 

