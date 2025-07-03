function B = set_up_matrix(G, mesh)
% B = set_up_matrix(G, mesh) 

% Load parameters
N = G.N;
Cm = G.Cm;
dt = G.dt;
Me = G.Me;
Mm = G.Mm;
chi = reshape(G.chi, G.N, 1);

% Set up sub matrices
I = speye(N,N);
Bm = set_up_sub_matrix(G, mesh, Mm);
Be = set_up_sub_matrix(G, mesh, Me);

% Combine matrices
B_top = [I+(dt/Cm)*Bm./chi, (dt/Cm)*Bm./chi];
B_top = B_top(G.with_cell,:);
B_top = B_top(:, [G.with_cell; G.N+(1:G.N)']);
B_bottom = [Bm, Bm+Be];

% Apply Dirichlet bc
for i=mesh.boundary
    last_row = zeros(1, 2*G.N);
    last_row(G.N+i) = 1;
    B_bottom(i,:) = last_row;
end
              
B_bottom = B_bottom(:, [G.with_cell; G.N+(1:G.N)']);                
B = [B_top; B_bottom];
B = sparse(B);

end
