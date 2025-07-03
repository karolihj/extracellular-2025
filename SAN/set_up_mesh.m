function mesh = set_up_mesh(G)
%mesh = set_up_mesh(G) Set up vectors containing the indices of the different 
%types of nodes

% Read geometry
Nx = G.Nx;
Ny = G.Ny;
Nz = G.Nz;
N = G.N;


% Set up indices for outer boundary
% Set up indices for outer boundary (corners)
mesh.e_lsw = 1;
mesh.e_lse = Nx;
mesh.e_lnw = (Ny-1)*Nx + 1;
mesh.e_lne = Nx*Ny;
mesh.e_hsw = (Nz-1)*Ny*Nx+1;
mesh.e_hse = (Nz-1)*Ny*Nx+Nx;
mesh.e_hnw = ((Nz-1)*Ny+Ny-1)*Nx + 1;
mesh.e_hne = Nx*Ny*Nz;

% Set up indices for outer boundary (lines)
mesh.e_lw = ((2:Ny-1)-1)*Nx + 1;
mesh.e_le = ((2:Ny-1)-1)*Nx + Nx;
mesh.e_ls = 2:Nx-1;
mesh.e_ln = (Ny-1)*Nx + (2:Nx-1);
mesh.e_hw = ((Nz-1)*Ny+(2:Ny-1)-1)*Nx + 1;
mesh.e_he = ((Nz-1)*Ny+(2:Ny-1)-1)*Nx + Nx;
mesh.e_hs = (Nz-1)*Ny*Nx+(2:Nx-1);
mesh.e_hn = ((Nz-1)*Ny+Ny-1)*Nx + (2:Nx-1);
mesh.e_sw = (1:Nz-2)*Ny*Nx+1;
mesh.e_se = (1:Nz-2)*Ny*Nx+Nx;
mesh.e_nw = ((1:Nz-2)*Ny+Ny-1)*Nx+1;
mesh.e_ne = ((1:Nz-2)*Ny+Ny-1)*Nx+Nx;

% Set up indices for outer boundary (sides)
[x, y] = meshgrid(2:Nx-1, 2:Ny-1);
mesh.e_l = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, ones((Ny-2),(Nx-2))), (Nx-2)*(Ny-2), 1))';
mesh.e_h = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, Nz*ones((Ny-2),(Nx-2))), (Nx-2)*(Ny-2), 1))';
[x, z] = meshgrid(2:Nx-1, 2:Nz-1);
mesh.e_s = sort(reshape(sub2ind([Nx,Ny,Nz], x, ones((Nz-2), (Nx-2)), z), (Nx-2)*(Nz-2), 1))';
mesh.e_n = sort(reshape(sub2ind([Nx,Ny,Nz], x, Ny*ones((Nz-2), (Nx-2)), z), (Nx-2)*(Nz-2), 1))';
[y, z] = meshgrid(2:Ny-1, 2:Nz-1);
mesh.e_w = sort(reshape(sub2ind([Nx,Ny,Nz], ones((Nz-2), (Ny-2)), y, z), (Ny-2)*(Nz-2), 1))';
mesh.e_e =  sort(reshape(sub2ind([Nx,Ny,Nz], Nx*ones((Nz-2), (Ny-2)), y, z), (Ny-2)*(Nz-2), 1))';

% Set up indices for the remaining extracellular cells
indices = ones(1, N);
indices([mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, mesh.e_hsw, mesh.e_hse, ...
            mesh.e_hnw, mesh.e_hne, mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, ...
            mesh.e_lw, mesh.e_le, mesh.e_ls, mesh.e_ln, mesh.e_ne, mesh.e_sw, ...
            mesh.e_se, mesh.e_nw, mesh.e_w, mesh.e_e, mesh.e_s, mesh.e_n, mesh.e_h, ...
            mesh.e_l]) = 0;
mesh.e = round(find(indices));

mesh.boundary = [mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, mesh.e_hsw, mesh.e_hse, ...
            mesh.e_hnw, mesh.e_hne, mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, ...
            mesh.e_lw, mesh.e_le, mesh.e_ls, mesh.e_ln, mesh.e_ne, mesh.e_sw, ...
            mesh.e_se, mesh.e_nw, mesh.e_w, mesh.e_e, mesh.e_s, mesh.e_n, mesh.e_h, ...
            mesh.e_l];


end