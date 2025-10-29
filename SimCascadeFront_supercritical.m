%% script to simulate the interaction of front interfaces in a system of two coupled FKPP equations
% the pde dynamics is simulated using a superposition of two supercritical fronts as initial condition
% the front speeds are then extracted from the simulation results
% the simulation uses MATLAB's pdepe solver

clear; close all; clc;

genVideo = 1; % set to 1 to generate video of the simulation

d = 4;
r = 2;
a1 = 0.75;
a2 = 0.75;

xmin = -500;
xmax = 1000;

%% calculate single front profiles

e4 = [(r+a1)/(r-a1*a2); 0; r*(1+a2)/(r-a1*a2); 0]; % left equilibrium
e1 = [1; 0; 0; 0]; % middle equilibrium
e0 = [0; 0; 0; 0]; % trivial equilibrium

% faster front
c2 = 2*sqrt(r*d) + 10; % set speed slightly
disp(['Set speed of faster front to c_2 = ',num2str(c2)]);

vecfield = @(~,u) [u(2);...
                   -(1/d)*(c2*u(2) + r*u(1).*(1-u(1)))];

L1 = [0,1;...
      r/d, -c2/d];

[eigvecs, eigvals] = eig(L1);
[sorted_eigenvals, sortidx] = sort(diag(eigvals));
sorted_eigenvecs = real(eigvecs(:,sortidx));

u_init = [1;0] + 1e-5*sorted_eigenvecs(:,2)*(-sorted_eigenvecs(1,2));
[xi2, sol_front2] = ode15s(vecfield, [0,xmax], u_init);
sol_front2 = sol_front2';

L0 = [0,1;...
      -r/d, -c2/d];
[eigvecs0, eigvals0] = eig(L0);
[sorted_eigenvals0, sortidx0] = sort(diag(eigvals0));
sorted_eigenvecs0 = real(eigvecs0(:,sortidx0));

xi_replace2_idx = find(sol_front2(1,:)<1e-4,1);
xi_refined = linspace(xi2(xi_replace2_idx), xi2(end), 5000);
sol_front2(:,xi_replace2_idx:(xi_replace2_idx+5000-1)) = sol_front2(:,xi_replace2_idx).*exp(-sorted_eigenvals0(2)*xi2(xi_replace2_idx)).*exp(sorted_eigenvals0(2)*xi_refined).*sorted_eigenvecs0(:,2);
xi2 = [xi2(1:xi_replace2_idx-1); xi_refined'];
sol_front2 = [sol_front2; zeros(2,length(xi2))]; % add zeros for u2 components

plot(xi2, sol_front2(1,:), 'Linewidth', 2);

% [~, xi2tilde, sol_front2tilde] = SolveBVP_coupled_front(c2, d, r, a1, a2, [xmin,xmax], 1500, e1, e0);

% slower front
c1 = 0.5*c2; % set speed slightly below c2
disp(['Set speed of slower front to c_1 = ',num2str(c1)]);

L1full = [0,1,0,0;...
          r/d, -c1/d, -a1, 0;...
          0,0,0,1;...
          0,0,-1-a2, -c1];
[eigvecs1, eigvals1] = eig(L1full);
[sorted_eigenvals1, sortidx1] = sort(diag(eigvals1));
sorted_eigenvecs1 = real(eigvecs1(:,sortidx1));

[~, xi1, sol_front1] = SolveBVP_coupled_front(c1, d, r, a1, a2, [xmin,0], 5000, e4, e1);

xi_replace1_idx = find(sol_front1(1,:)-1<1e-1,1);
sol_front1(:,xi_replace1_idx:end) = e1 + (sol_front1(:,xi_replace1_idx)-e1.*ones(size(sol_front1(:,xi_replace1_idx)))).*exp(-sorted_eigenvals1(3)*xi1(xi_replace1_idx)).*exp(sorted_eigenvals1(3)*xi1(xi_replace1_idx:end)).*sorted_eigenvecs1(:,3);

xgrid = linspace(xmin,xmax,2000);
u_init = frontSuperposition(xgrid,0,xi1,sol_front1,xi2,sol_front2,c2,c1);
figure(4);
plot(xgrid,u_init(1,:),'Linewidth',2);
hold on;
plot(xgrid,u_init(2,:),'Linewidth',2);
hold off;

%% pde dynamics

m = 0;
tend = 30;
x = xmin:0.05:xmax;
t = 0:0.1:tend;

tic;
sol = pdepe(m,@pdex4pde,@(x)pdex4ic(x,xgrid,u_init),@(xl,ul,xr,ur,t,er,el)pdex4bc(xl,ul,xr,ur,t,e4,e0),x,t);
% sol = pdepe(m,@(x,t,u,DuDx)pdex4pdeLinear(x,t,u,DuDx,xi1,sol_front1,xi2,sol_front2,c2,c1),@(x)pdex4ic(x,xgrid,u_init),@(xl,ul,xr,ur,t,er,el)pdex4bc(xl,ul,xr,ur,t,e4,e0),x,t);
runTime = toc;
disp(['Run time: ',num2str(runTime)]);
u1 = sol(:,:,1);
u2 = sol(:,:,2);

%% calculate selected front speeds

[c1sel,c2sel,f1init,f2init,t0] = frontSpeed(u1,u2,t,x);

disp(['Speed of leading front: c_1 = ',num2str(c1sel)]);
disp(['Speed of secondary front: c_2 = ',num2str(c2sel)]);

%% 
if genVideo
    v = VideoWriter('front-cascade-supercrit-init.mp4','MPEG-4');
    open(v);
    fig = figure;
    h = figure(5);
    for i = 1:length(t)
        % X = linspace(-20,20,500);
        plot(x,u1(i,:),'Linewidth',2);
        axis([xmin+10,xmax,-5,5]);
        % axis off;
        hold on;
        plot(x,u2(i,:),'Linewidth',2);
        % plot(frontpos1(i),1.5,'x');
        % axis off;
        caption = sprintf('t=%f',t(i));
        title(caption, 'FontSize', 20);
        hold off;
        writeVideo(v,getframe(h));
    end
    close(v);
end


%% plot solution in space-time plot with top-down view

[X,T] = meshgrid(x,t);
figure(2);
hold on
s = pcolor(X,T,u1);
% view(2);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
xlim([xmin+10,xmax])
ylim([0,tend])
pbaspect([100/80 1 1])

set(gca, 'XTick', [], 'YTick', [])
ax = gca;                 % Get current axes
ax.FontSize = 36;         % Change font size
ax.LineWidth = 2;

annotation('arrow', [0.148,0.9], [0.133 0.133], 'LineWidth', 2) % x-axis arrow
annotation('arrow', [0.148 0.148], [0.13 0.94], 'LineWidth', 2) % y-axis arrow

% plot speed of first front
xidx = find(x>f1init,1);
plot(x(xidx:end),(x(xidx:end)-f1init)/c1sel + t0,'Color','red','LineWidth',5)

% plot speed of second front
xidx = find(x>f2init,1);
plot(x(xidx:end),(x(xidx:end)-f2init)/c2sel + t0,'Color','green','LineWidth',5)
hold off
exportgraphics(gca,'front-cascade-supercrit-init.jpg','Resolution',600)

%% plot solution in space-time plot with side view
 
figure(3);
hold on;
for ii = 1:10:numel(t)
    plot3(x,t(ii)*ones(size(x)),u1(ii,:),'k');
    xlim([xmin+10,xmax]);
end
view(15,20)
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
ax = gca;                 % Get current axes
ax.FontSize = 15;         % Change font size
% ax.LineWidth = 2;
pbaspect([100/30 1 1])
exportgraphics(gca,'front-cascade-supercrit-init.jpg','Resolution',600)


%% --------------------------------------------------------------------------
% define functions


function [c,f,s] = pdex4pde(x,t,u,DuDx) % sets up the pde for the simulation

% set parameters
d = 4;
r = 2;
a1 = 0.75;
a2 = 0.75;

% define pde
c = [1;1];
f = [d 0; 0 1] * DuDx;
s = [r 0; 0 1]*u.*(ones(size(u))-u) + ([a1 0; 0 a2]*u).*([0 1;1 0]*u);
end

% --------------------------------------------------------------------------

function [c,f,s] = pdex4pdeLinear(x,t,u,DuDx,xi1,front1,xi2,front2,c2,c1) % sets up the pde for the simulation of linear dynamics

% set parameters
d = 4;
r = 2;
a1 = 0.75;
a2 = 0.75;

uSuperpos = frontSuperposition(x,t,xi1,front1,xi2,front2,c2,c1);

% define pde
c = [1;1];
f = [d 0; 0 1] * DuDx;
s = [r 0; 0 1]*u.*(1-2*uSuperpos) + ([a1 0; 0 a2]*uSuperpos).*([0 1;1 0]*u) + ([a1 0; 0 a2]*u).*([0 1;1 0]*uSuperpos);
end

% --------------------------------------------------------------------------

function u0 = pdex4ic(x,xgrid,u_init) % sets initial profile

u0 = [interp1(xgrid,u_init(1,:),x,'nearest','extrap');...
      interp1(xgrid,u_init(2,:),x,'nearest','extrap')];

% u0 = 1e-5*(x>100).*(x<120)*[1;1]; % small perturbation initial condition

end
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t,er,el) % implements dirichlet boundary conditions
uf1 = el(1);
uf2 = el(2);

pl = ul-[uf1;uf2];
ql = [0; 0];
pr = ur;
qr = [0; 0];

end
% --------------------------------------------------------------------------


function [c1,c2,front1init,front2init,t0] = frontSpeed(u1,u2,t,x)
    startidx = floor(numel(t)/3); % start measurement only at half-way point
    t0 = t(startidx);

    frontpos2 = zeros(numel(t),1);
    frontpos1 = zeros(numel(t),1);
    
    for ii = startidx:numel(t)
        temp = u2(ii,:)>0.5;
        if sum(temp) > 0.5
            frontpos2(ii) = x(find(temp,1,'last'));
        end
    
        temp = u1(ii,:)>0.5;
        if sum(temp) > 0.5
            frontpos1(ii) = x(find(temp,1,'last'));
        end
    end

    front1init = frontpos1(startidx);
    front2init = frontpos2(startidx);
    
    speed2 = zeros(numel(t),1);
    speed1 = zeros(numel(t),1);


    for ii = startidx:numel(t)-1
        speed2(ii) = (frontpos2(ii+1)-frontpos2(ii))/(t(ii+1)-t(ii));
        speed1(ii) = (frontpos1(ii+1)-frontpos1(ii))/(t(ii+1)-t(ii));
    end
   
    c1 = sum(speed1(startidx:numel(t)-1))/numel(startidx:numel(t)-1);
    c2 = sum(speed2(startidx:numel(t)-1))/numel(startidx:numel(t)-1);
end
% --------------------------------------------------------------------------


function [u] = frontSuperposition(x,t,xi1,sol_front1,xi2,sol_front2,c2,c1)
    u = zeros(2,length(x));
    xcut = - (c2+c1)/2*t; % initial cut position at t
    u(1,:) = (x>xcut).*interp1(xi2, sol_front2(1,:), x-c2*t, 'nearest', 'extrap') + (x < xcut).*interp1(xi1, sol_front1(1,:), x-c1*t, 'nearest', 'extrap');
    u(2,:) = (x>xcut).*interp1(xi2, sol_front2(3,:), x-c2*t, 'nearest', 'extrap') + (x < xcut).*interp1(xi1, sol_front1(3,:), x-c1*t, 'nearest', 'extrap');
    % u(1,:) = interp1(xi2, sol_front2(1,:), x, 'nearest', 'extrap');
    % u(2,:) = interp1(xi2, sol_front2(3,:), x, 'nearest', 'extrap');
end
% --------------------------------------------------------------------------
