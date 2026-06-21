function plot_K(x, y, z)

tri = delaunay(x, y);
trisurf(tri, x, y, z);
shading interp; view(3); grid on; colorbar
xlabel x; ylabel y; zlabel z

end

