function plot_rho(mesh_divide_elements, mesh_divide_nodes, rho)
element_num = length(mesh_divide_elements);

for e = 1: element_num
    e_node = mesh_divide_elements(e);
    x_e = [mesh_divide_nodes(e_node(1), 1);...
           mesh_divide_nodes(e_node(2), 1);...
           mesh_divide_nodes(e_node(3), 1)];
    y_e = [mesh_divide_nodes(e_node(1), 2);...
           mesh_divide_nodes(e_node(2), 2);...
           mesh_divide_nodes(e_node(3), 2)];
    rho_e = [rho(e_node(1)); rho(e_node(2)); rho(e_node(3))];
    
    
    
end

