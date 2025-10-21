%% script to simulate the interaction of front interfaces in a system of two coupled FKPP equations
% the pde dynamics is simulated from compactly supported step initial data
% around x = 0 on a sufficiently large spacial interval such that boundary
% effects can be neglected.
% The solution is numerically calculated using pdepe.

genVideo = 1;


m = 0;
xmin = -800;
xmax = 800;
tend = 120;
x = xmin:0.1:xmax;
t = 0:0.1:tend;

d = 4;
r = 2;
a1 = 0.75;
a2 = 0.75;

tic;
sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t);
runTime = toc;
disp(['Run time: ',num2str(runTime)]);
u1 = sol(:,:,1);
u2 = sol(:,:,2);
%% calculate and find optimal weights in assumption 3

[c1,c2,f1init,f2init,t0] = frontSpeed(u1,u2,t,x);

disp(['Speed of leading front: c_1 = ',num2str(c1)]);
disp(['Speed of secondary front: c_2 = ',num2str(c2)]);
  
% define function handle
kappa = linspace(0,10,1000);
hold on
plot(kappa,kappa.^2-c2*kappa+(1+a2))
plot(kappa,d*kappa.^2-c2*kappa-r)
ylim([-5,1]);
hold off

pause()

plot(kappa,d*kappa.^2-c1*kappa+r)
ylim([-5,1]);

pause()

[kappa1,kappa2] = meshgrid(linspace(0,2,200),linspace(0,2,200));
surf(kappa1,kappa2,d*kappa1.^2-c1*kappa1+r+kappa2*(c1-c2));
view(2);

%% 
if genVideo
    v = VideoWriter('CompactInit.mp4','MPEG-4');
    open(v);
    fig = figure;
    h = figure(1);
    for i = 1:length(t)
        % X = linspace(-20,20,500);
        plot(x,u1(i,:),'Linewidth',2);
        axis([xmin,xmax,-5,5]);
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
xlim([xmin,xmax])
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
plot(x(xidx:end),(x(xidx:end)-f1init)/c1 + t0,'Color','red','LineWidth',5)

% plot speed of second front
xidx = find(x>f2init,1);
plot(x(xidx:end),(x(xidx:end)-f2init)/c2+ t0,'Color','green','LineWidth',5)
hold off
exportgraphics(gca,'front-cascade-compact-init-top-view.jpg','Resolution',600)

%% plot solution in space-time plot with side view
 
figure(3);
hold on;
for ii = 1:10:numel(t)
    plot3(x,t(ii)*ones(size(x)),u1(ii,:),'k');
    xlim([0,xmax]);
end
view(15,20)
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
ax = gca;                 % Get current axes
ax.FontSize = 20;         % Change font size
% ax.LineWidth = 2;
pbaspect([100/30 1 1])
exportgraphics(gca,'front-cascade-compact-init-side-view.jpg','Resolution',600)


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

function u0 = pdex4ic(x) % sets initial profile
 
% eps = 1e-2;
% u0 = eps*(x<5)*(x>-5)*[0;1]+(x<40)*(x>-40)*[1;0];
u0 = 0.1*(x<10)*(x>-10)*[1;1];

end
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t) % implements zero dirichlet boundary conditions
uf1 = 0;
uf2 = 0;

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