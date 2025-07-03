function B = set_up_matrix(G, mesh)
% B = set_up_matrix(G, mesh) 

% Load parameters
N = G.N;
Cm = G.Cm;
dt = G.dt;
Me = G.Me;
Mm = G.Mm;
Am = G.Am(:);

% Set up sub matrices
I = speye(N,N);
Bm = set_up_sub_matrix(G, mesh, Mm);
Be = set_up_sub_matrix(G, mesh, Me);

% Combine matrices
B_top = [I+(dt/Cm)*Bm./Am, (dt/Cm)*Bm./Am]; 
B_bottom = [Bm, Bm+Be];

% Apply Dirichlet bc
for i=[mesh.e_s, mesh.e_n, mesh.e_w, mesh.e_e, mesh.e_sw, mesh.e_se, mesh.e_ne, mesh.e_nw]
	last_row = zeros(1, 2*N);
	last_row(G.N+i) = 1;
    B_bottom(i,:) = last_row;
end

B_top = B_top(G.with_myocyte,:);
B_top = B_top(:, [G.with_myocyte; G.N+(1:G.N)']);
B_bottom = B_bottom(:, [G.with_myocyte; G.N+(1:G.N)']);
B = [B_top; B_bottom];

B = sparse(B);
end
