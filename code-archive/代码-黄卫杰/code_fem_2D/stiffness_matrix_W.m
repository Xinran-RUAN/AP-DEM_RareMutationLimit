function [A, b] = stiffness_matrix_W(K1, K2_x, K2_y, Con, mesh_divide_elements,...
                            mesh_divide_nodes, xi, eta, weight,...
                            D, W, du_m_dw, K, rho, H, eps, dt, d_theta, N_theta)
element_num = length(mesh_divide_elements);
node_num = length(mesh_divide_nodes);

A = zeros(node_num * N_theta);
b = zeros(node_num * N_theta, 1);

for kk = 1: N_theta
    D_k = D(kk); 
    W_k = W(:, kk);
    H_k = H(kk);
    
    du_m_dw_k = du_m_dw(:, kk);
    
    for e = 1: element_num
        
        e_node = mesh_divide_elements(e, :);
        
        W_local = [W_k(e_node(1)); W_k(e_node(2)); W_k(e_node(3))];
        du_m_dw_local = [du_m_dw_k(e_node(1));...
                         du_m_dw_k(e_node(2)); du_m_dw_k(e_node(3))];
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
  
        [A1, A2, A3, A4, b1, b2] = element_stiffness_matrix_W(K1, k2_x_e, k2_y_e,...
                               Con, x_local, y_local, xi, eta, weight, D_k, H_k,...
                               W_local, du_m_dw_local, K_local, rho_local,...
                                                            eps, dt, d_theta);
    
        nd1 = e_node + node_num * (kk - 1);
        nd2 = e_node + node_num * (kk - 1);
        A(nd1, nd2) = A(nd1, nd2) + A1 + A2 - 2 * A3 + A4;
        
        if kk == N_theta
            nd3 = nd2 - node_num * (N_theta - 1);
        else
            nd3 = nd2 + node_num;
        end
        A(nd1, nd3) = A(nd1, nd3) + A3;
        
        if kk == 1
            nd4 = nd2 + node_num * (N_theta - 1);
        else
            nd4 = nd2 - node_num;
        end
        A(nd1, nd4) = A(nd1, nd4) + A3;
        
        b(nd1) = b(nd1) + b1 + b2;
    end
end 
  
end

