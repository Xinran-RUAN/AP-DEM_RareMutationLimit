function ntheta = int_theta_n(mesh_divide_elements, mesh_divide_nodes, n, N_theta)

element_num = length(mesh_divide_elements);
ntheta = zeros(1, N_theta);

for kk = 1: N_theta
    n_k = n(:, kk);
    for  e = 1:element_num
        e_node = mesh_divide_elements(e, :);
        x_local = [mesh_divide_nodes(e_node(1), 1);...
            mesh_divide_nodes(e_node(2), 1);...
            mesh_divide_nodes(e_node(3), 1)];
        y_local = [mesh_divide_nodes(e_node(1), 2);...
            mesh_divide_nodes(e_node(2), 2);...
            mesh_divide_nodes(e_node(3), 2)];
        
        n_local = [n_k(e_node(1)); n_k(e_node(2)); n_k(e_node(3))];
        Jacobi = (x_local(2) - x_local(1)) * (y_local(3) - y_local(1)) - (x_local(3) - x_local(1)) * (y_local(2) - y_local(1));
        
        ntheta(kk) = ntheta(kk) +...
            + 1/6 * (n_local(1) + n_local(2) + n_local(3)) * Jacobi;
    end
end
end

