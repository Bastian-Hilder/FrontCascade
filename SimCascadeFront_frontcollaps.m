%% script to simulate the interaction of front interfaces in a system of two coupled FKPP equations
% the pde dynamics is simulated from compactly supported step initial data
% around x = 0 on a sufficiently large spacial interval such that boundary
% effects can be neglected.
% The solution is numerically calculated using pdepe.

clear all;
clc;
close all;

genVideo = 0;


m = 0;
xmin = -120;
xmax = 120;
tend = 40;
x = xmin:0.1:xmax;
t = 0:0.1:tend;

tic;
sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t);
runTime = toc;
disp(['Run time: ',num2str(runTime)]);
u1 = sol(:,:,1);
u2 = sol(:,:,2);
%% calculate and plot front speeds

[f1pos,f2pos,t0idx] = frontPos(u1,u2,t,x);


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
xlabel("x")
ylabel("t")
ylim([0,tend])
xlim([xmin,xmax])
pbaspect([100/30 1 1])

plot(f1pos,t(t0idx:end),'Color','red','LineWidth',2)

plot(f2pos,t(t0idx:end),'Color','green','LineWidth',2)
hold off
exportgraphics(gca,'front-cascade-front-collaps-top-view.jpg','Resolution',600)

%% plot solution in space-time plot with side view

figure(3);
hold on;
for ii = 1:10:numel(t)
    plot3(x,t(ii)*ones(size(x)),u1(ii,:),'k');
    xlim([0,xmax]);
end
view(15,20)
xlabel("x")
ylabel("t")
pbaspect([100/30 1 1])
exportgraphics(gca,'front-cascade-front-collaps-side-view.jpg','Resolution',600)


%% --------------------------------------------------------------------------
% define functions


function [c,f,s] = pdex4pde(x,t,u,DuDx) % sets up the pde for the simulation

% set parameters
d = 0.2;
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

d = 0.2;
r = 2;
a1 = 0.75;
a2 = 0.75;

e1 = [(r+a1)/(r-a1*a2);(r*(1+a2))/(r-a1*a2)];

u0 = (x<50)*(x>-50)*[1;0] + (x<10)*(x>-10)*(e1-[1;0]);

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


function [frontpos1,frontpos2,t0idx] = frontPos(u1,u2,t,x)
    startidx = floor(numel(t)/4); % start measurement only at half-way point
    t0idx = startidx;

    frontpos2 = zeros(numel(t),1);
    frontpos1 = zeros(numel(t),1);
    
    for ii = startidx:numel(t)
        temp = u1(ii,:)>1.5;
        if sum(temp) > 0.5
            frontpos2(ii) = x(find(temp,1,'last'));
        end
    
        temp = u1(ii,:)>0.5;
        if sum(temp) > 0.5
            frontpos1(ii) = x(find(temp,1,'last'));
        end
    end
    
    speed2 = zeros(numel(t),1);
    speed1 = zeros(numel(t),1);


    for ii = startidx:numel(t)-1
        speed2(ii) = (frontpos2(ii+1)-frontpos2(ii))/(t(ii+1)-t(ii));
        speed1(ii) = (frontpos1(ii+1)-frontpos1(ii))/(t(ii+1)-t(ii));
    end

    frontpos1 = frontpos1(startidx:end);
    frontpos2 = frontpos2(startidx:end);
end