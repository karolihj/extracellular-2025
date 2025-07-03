function [P] = vary_parameters(P, G)
%[P] = vary_parameters(P, G)

[x, y, z] = meshgrid(1:G.Nx, 1:G.Ny, 1:G.Nz);
x = permute(x,[2 1 3]);
y = permute(y,[2 1 3]);
x = x(:);
y = y(:);

x = x(G.with_cell);
y = y(G.with_cell);

x_th = round(G.Nx/2);
y_th = round(G.Ny/2);

change_idx = [81, 12, 84, 88, 89, 86, 87, 62, 66, 73];

for n=1:G.Nm
    if (x(n) <= x_th) && (y(n) > y_th)
        P(change_idx, n) = P(change_idx, n)*(1+G.vary_percentage/2);
    elseif (x(n) > x_th) && (y(n) <= y_th)
        P(change_idx, n) = P(change_idx, n)*(1-G.vary_percentage/2);
    end
end


end