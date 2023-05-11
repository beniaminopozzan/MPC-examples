function [u_opt, U, X] = LinearMPC(plant, x0, U0, ref, MPC)


%----------OPC2QP----------
A = plant.A;
B = plant.B;

% n state size, m input size
[nx,nu] = size(B);

N     = MPC.N;                    
Q     = MPC.Q;    
R     = MPC.R;
P     = MPC.P;
u_min = MPC.u_min;
u_max = MPC.u_max;
x_max = MPC.x_max;

% A_bar generation
A_bar = zeros(nx*N,nx);           %initialize A_bar
A_bar(1:nx,:) = A;
A_power = A;                    %auxiliary variable for generation of A_bar 
for i = 2:N
    A_power = A_power*A;
    A_bar((i-1)*nx+1:i*nx,:) = A_power; 
end

% B_bar generation
B_bar=zeros(nx*N,nu*N);           %initialize B_bar
B_bar_row = [B zeros(nx,(N-1)*nu)]; 
B_bar(1:nx,:) = B_bar_row; 
for i = 2:N
    B_bar_row = [A*B_bar_row(1:nx,1:nu),B_bar_row(:,1:(N-1)*nu)]; 
    B_bar((i-1)*nx+1:i*nx,:) =  B_bar_row;
end

% R_bar generation
R_bar = kron(eye(N),R);

% Q_bar generationt
Q_bar = blkdiag(kron(eye(N-1),Q),P);

X_ref = repmat(ref,N,1);

% Quadprog matrices generation 
H = R_bar + B_bar'*Q_bar*B_bar;
H =(H+H')/2;
f = (A_bar*x0 - X_ref)'*Q_bar*B_bar;

ub = repmat(u_max,N,1);
lb = repmat(u_min,N,1);

% state const
A_con = B_bar;
b_con = repmat(x_max,N,1) - A_bar*x0;

%----------Quadratic programming solution----------
options = optimoptions('quadprog','Display','off','Algorithm', 'interior-point-convex');

%U0 = zeros(nu*N,1);
[U, ~, exitflag] = quadprog(H,f,A_con,b_con,[],[],lb,ub,U0,options);

if( exitflag ~= 1 )
    error("failed to solve the QP problem")
end

u_opt = U(1:nu);

X = reshape(A_bar*x0+B_bar*U,nx,N);
U = reshape(U,nu,N);