function mesh_triangle(iteration_max, h, xa, xb, yc, yd)
if ( nargin < 1)
    iteration_max = 200;
    fprintf ( 1, '\n');
    fprintf ( 1, ' mesh_divide-Note:\n');
    fprintf ( 1, ' No value of ITERATION_MAX was supplied.\n');
    fprintf ( 1, ' The default value ITERATION_MAX = %d will be used.\n', ...
        iteration_max);
end

if ( nargin < 2)
    h = 0.005;
    fprintf ( 1, '\n');
    fprintf ( 1, ' mesh_divide-Note:\n');
    fprintf ( 1, ' No value of H was supplied.\n');
    fprintf ( 1, ' The default value H = %f will be used.\n', h);
end

if ( nargin < 3)
    fh = @mesh_divide_pie_slice_fh;
    fprintf ( 1, '\n');
    fprintf ( 1, ' mesh_divide-Note:\n');
    fprintf ( 1, ' No value of FH was supplied.\n');
    fprintf ( 1, '  The variable density function FHHOLEY = %f will be used.\n', h );
end

% put the random number generator into a fixed initial state.

rand('state', 111);

% set the rendering method for the current figure to Z-buffering

set( gcf, 'rend', 'z');

fprintf ( 1, '\n');
fprintf ( 1, ' Problem holey pie slice:\n');
fprintf ( 1, ' The holey pie slice, h = %f\n', h)

fd = @(p) drectangle (p,xa,xb,yc,yd);
% fh = @(p) pro_holey_fh(p,cc1,cc2,zeta);
box = [xa-1, yc-1; xb+1, yd+1];
xf = [xa: 0.5: xb]; 
yf = [yc: 0.5: yd];
fixed = zeros(length(xf) * length(yf), 2);
for kk = 1: length(xf)
    for kkk = 1: length(yf)
        num = (kkk - 1) * length(xf) + kk;
        fixed(num, :) = [xf(kk), yf(kkk)];
    end
end

% fixed = [];
[p,t] = distmesh2d (fd, @huniform, h, box, fixed);

% post_2d ( p, t, fh)

% write a post script image of the triangulation.

[ node_num, junk] = size (p);
[ tri_num, junk] = size (t);
p = p';
t = t';
node_show = 0;
triangle_show = 1;

triangulation_order3_plot( 'mesh_divide_mesh.eps', node_num, p, tri_num, ...
    t, node_show, triangle_show);

% write a text file containing the nodes.

r8mat_write ( 'mesh_divide_nodes.txt', 2, node_num, p);

% write a text file containint the triangles.

i4mat_write( 'mesh_divide_elements.txt', 3, tri_num, t);

triangulation_boundary_nodes('mesh_divide');

triangulation_l2q('mesh_divide');

triangulation_boundary_nodes('mesh_divide_l2q');

return
end

