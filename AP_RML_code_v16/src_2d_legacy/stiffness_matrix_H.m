function [A, B] = stiffness_matrix_H(mesh_divide_elements, mesh_divide_nodes, D, K, rho, xi, eta, weight, K1, K2_x, K2_y)
element_num = length(mesh_divide_elements);
node_num = length(mesh_divide_nodes);

A = zeros(node_num);
B = zeros(node_num);

for e = 1: element_num
 
    e_node = mesh_divide_elements(e, :);
    
    K_local = [K(e_node(1)); K(e_node(2)); K(e_node(3))];
    rho_local = [rho(e_node(1)); rho(e_node(2)); rho(e_node(3))];
    k2_x_e = K2_x(:, :, :, e);
    k2_y_e = K2_y(:, :, :, e);
    
    x_local = [mesh_divide_nodes(e_node(1), 1);...
               mesh_divide_nodes(e_node(2), 1);...
               mesh_divide_nodes(e_node(3), 1)];
    y_local = [mesh_divide_nodes(e_node(1), 2);...
               mesh_divide_nodes(e_node(2), 2);...
               mesh_divide_nodes(e_node(3), 2)];
    [A_local, B_local] = element_stiffness_matrix_H(K1, k2_x_e, k2_y_e,...
                         x_local, y_local, K_local, rho_local, xi, eta, weight, D);
   
    nd1 = e_node;
    nd2 = e_node;
    A(nd1, nd2) = A(nd1, nd2) + A_local;
    B(nd1, nd2) = B(nd1, nd2) + B_local;
end
end

