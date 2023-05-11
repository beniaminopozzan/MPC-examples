%% Model Parameters

% model in state space of a 2D double integrator
A = [
    0 0 1 0 
    0 0 0 1
    0 0 0 0
    0 0 0 0
    ];
B = [
    0 0
    0 0
    1 0
    0 1
    ];
C = diag([ ...
    1 1 0 0 ...
    ]);
D = [
    0 0
    0 0
    0 0
    0 0
    ];
plant = ss(A, B, C, D);
[nx, nu] = size(B);

%% MPC and Simulation Parameters

% sampling time
Ts = 0.01;

% MPC parameters
MPC.N     = 50;                 
MPC.Q     = diag([50 50 1 1]);
MPC.R     = 1*diag([1 1]);
MPC.P     = MPC.Q;
MPC.u_min = -1*[1 1]'*inf;
MPC.u_max = 1*[1 1]'*inf;
MPC.x_max = [1.2 2.2 2 2]'*inf;

% Simulation Parameters

x0 = [
    0
    0
    0
    0
    ];
ref = [
    1
    2
    0
    0
    ];
sim_t = 10;


%% MPC Simulation

% plant discretization
plant_d = c2d(plant, Ts, 'zoh');

% variable allocation
t = (0:sim_t/Ts)*Ts;
T = length(t);
x = zeros(nx,T);
u = zeros(nu,T-1);

x(:,1) = x0;

t_ext = (0:sim_t/Ts+MPC.N)*Ts;
T_ext = length(t_ext);
u_vec = repmat([0 0]',1,MPC.N);
x_vec = repmat(x0,1,MPC.N);

figure(1)

set(gcf,'DefaultLineLineWidth',1.5)
set(gcf,'DefaultAnimatedLineLineWidth',1.5)
set(gcf,'DefaultAxesLineWidth',2)
set(gcf,'DefaultAxesFontSize',15)

set(gcf,'defaultAxesTickLabelInterpreter','latex');
set(gcf,'defaultLegendInterpreter','latex');
set(gcf,'defaulttextinterpreter','latex');
set(gcf,'defaultGraphplotInterpreter','latex');


tl = tiledlayout(3,2);
title(tl,"trajectories", 'interpreter','latex');
ax(1) = nexttile(1);
ylabel("$p(t)$, [m]")
hold on
plot([0 t_ext(end)],(ref(1:2)*[1 1])',':')
ax(1).ColorOrderIndex = 1;
plot([0 t_ext(end)],(MPC.x_max(1:2)*[1 1])','-.')
ax(1).ColorOrderIndex = 1;
pred(1:2) = plot(t_ext(2:MPC.N+1),x_vec(1:2,:),'--');
x_and_u_lines(1:2) = arrayfun(@(x) animatedline(0,x),x0(1:2));
ylim([0 4]);

ax(2) = nexttile(3);
ylabel("$v(t)$, [m/s]")
hold on
plot([0 t_ext(end)],(ref(3:4)*[1 1])',':')
ax(2).ColorOrderIndex = 1;
plot([0 t_ext(end)],(MPC.x_max(3:4)*[1 1])','-.')
ax(2).ColorOrderIndex = 1;
pred(3:4) = plot(t_ext(2:MPC.N+1),x_vec(3:4,:),'--');
x_and_u_lines(3:4) = arrayfun(@(x) animatedline(0,x),x0(3:4));
ylim([-3 8]);

ax(3) = nexttile(5);
ylabel("$u(t)$, [m/s2]")
hold on
plot([0 t_ext(end)],([0 0]'*[1 1])',':')
ax(3).ColorOrderIndex = 1;
plot([0 t_ext(end)],(MPC.u_max*[1 1])','-.')
ax(3).ColorOrderIndex = 1;
plot([0 t_ext(end)],(MPC.u_min*[1 1])','-.')
ax(3).ColorOrderIndex = 1;
pred(5:6) = plot(t_ext(1:MPC.N),u_vec(1:2,:),'--');
x_and_u_lines(5:6) = arrayfun(@(x) animatedline(),[0 0]);
ylim([-8 8]);

ax(4) = nexttile(2,[3,1]);
hold on
axis equal
plot([1 1 0]*MPC.x_max(1), [0 1 1]*MPC.x_max(2),'-.k')
plot(x0(1),x0(2),'o')
plot(ref(1),ref(2),'*')
traj = animatedline(x0(1),x0(2),'color','k');
traj_pred = plot(x_vec(1,:),x_vec(2,:),':r');
xlim([0 4])
ylim([0 4])
xlabel("$x$, [m]");
ylabel("$y$, [m]");


set(ax,'XGrid','on')
set(ax,'YGrid','on')
set(ax(1:3),'XLim',[0 t_ext(end)])
arrayfun(@(axes) xlabel(axes,"$t$, [s]", 'interpreter','latex'), ax(1:3));

set(x_and_u_lines(1:2:end),'Color',"#0072BD")
set(x_and_u_lines(2:2:end),'Color',"#D95319")


for i=1:T-1
    [u_opt, u_vec, x_vec] = LinearMPC(plant_d, x(:,i), u_vec, ref, MPC);
    u(:,i) = u_opt;
    x(:,i+1) = plant_d.A*x(:,i) + plant_d.B*u_opt; % + chol(diag([0.001 0.001 0.001 0.0001]))*randn(4,1);
    u_vec = [u_vec(:,2:end) u_vec(:,end)];

    set(pred,'XData',t_ext(i-1+(2:MPC.N+1)));
    set(pred(1:4),{'YData'},num2cell(x_vec,2));
    set(pred(5:6),{'YData'},num2cell(u_vec,2));
    set(traj_pred,{'XData','YData'},num2cell(x_vec(1:2,:),2)')
    arrayfun(@(l,t,y) addpoints(l,t,y), x_and_u_lines', t([i+1 i+1 i+1 i+1 i i])', [x(:,i+1); u(:,i)])
    addpoints(traj, x(1,i+1), x(2,i+1))
    pause(Ts)
    if( i==1 )
        input("ready?")
    end
end