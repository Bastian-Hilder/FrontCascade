%% Boundary Value Problem Solver for Coupled Front System
% Solves the 4D ODE system for a traveling wave connecting two equilibria
% Uses MATLAB's bvp4c solver for boundary value problems

function [sol, xi_mesh, u_sol] = SolveBVP_coupled_front(c1, d, r, a1, a2, xi_span, num_points, el, er)
% Inputs:
%   c1         - wave speed
%   d          - diffusion coefficient ratio
%   r          - growth rate parameter
%   a1, a2     - coupling parameters
%   xi_span    - domain [xi_min, xi_max] for the traveling wave coordinate
%   num_points - number of mesh points for initial guess
%   el         - left equilibrium in the form [u1; u1'; u2; u2']
%   er         - right equilibrium in the form [u1; u1'; u2; u2']
%
% Outputs:
%   sol        - solution structure from bvp4c
%   xi_mesh    - mesh points where solution is computed
%   u_sol      - solution values [u1; u1'; u2; u2'] at mesh points

% Define the ODE system
odefun = @(xi, u) ode_system(xi, u, c1, d, r, a1, a2);

% Define boundary conditions
bcfun = @(ya, yb) boundary_conditions(ya, yb, el, er);

% Create initial mesh
xi_init = linspace(xi_span(1), xi_span(2), num_points);

% Create initial guess 
u_init = @(xi) [(el(1) - er(1))/2 * (tanh(-(xi+(xi_span(2)-xi_span(1))/2)) + 1) + er(1);...
                (el(2) - er(2))/2 * (tanh(-(xi+(xi_span(2)-xi_span(1))/2)) + 1) + er(2);...
                (el(3) - er(3))/2 * (tanh(-(xi+(xi_span(2)-xi_span(1))/2)) + 1) + er(3);...
                (el(4) - er(4))/2 * (tanh(-(xi+(xi_span(2)-xi_span(1))/2)) + 1) + er(4)];

% u1 = u_init(xi_init);

% plot(xi_init, u1(1,:), 'r'); % plot initial guess at first point


% Set up boundary value problem
solinit = bvpinit(xi_init, u_init);

% Solve the boundary value problem
options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Stats', 'on');
sol = bvp4c(odefun, bcfun, solinit, options);

% Extract solution
xi_mesh = sol.x;
u_sol = sol.y;

% Plot results
figure;
subplot(2,1,1)
plot(xi_mesh, u_sol(1,:), 'b-', 'LineWidth', 2);
hold on;
plot(xi_mesh, u_sol(3,:), 'r-', 'LineWidth', 2);
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$u_1, u_2$', 'Interpreter', 'latex', 'FontSize', 14);
legend('$u_1$', '$u_2$', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Traveling Wave Profile (c = %.4f)', c1));
grid on;

subplot(2,1,2)
plot(xi_mesh, u_sol(2,:), 'b--', 'LineWidth', 2);
hold on;
plot(xi_mesh, u_sol(4,:), 'r--', 'LineWidth', 2);
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$u_1'', u_2''$', 'Interpreter', 'latex', 'FontSize', 14);
legend('$u_1''$', '$u_2''$', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

% fprintf('\nBVP Solution Summary:\n');
% fprintf('  Wave speed c1 = %.6f\n', c1);
% fprintf('  Domain: [%.2f, %.2f]\n', xi_span(1), xi_span(2));
% fprintf('  Number of mesh points: %d\n', length(xi_mesh));
% fprintf('  Left equilibrium:  u1 = %.6f, u2 = %.6f\n', el(1), el(3));
% fprintf('  Right equilibrium: u1 = %.6f, u2 = %.6f\n', er(1), er(3));
% fprintf('  Solution at boundaries:\n');
% fprintf('    Left:  u1 = %.6e, u2 = %.6e\n', u_sol(1,1), u_sol(3,1));
% fprintf('    Right: u1 = %.6e, u2 = %.6e\n', u_sol(1,end), u_sol(3,end));

end

%% ODE System
function dudt = ode_system(~, u, c1, d, r, a1, a2)
% The 4D ODE system for the traveling wave
% u = [u1; u1'; u2; u2']
% dudt = [u1'; u1''; u2'; u2'']

dudt = zeros(4, 1);
dudt(1) = u(2);
dudt(2) = -(1/d) * (c1*u(2) + r*u(1)*(1-u(1)) + a1*u(1)*u(3));
dudt(3) = u(4);
dudt(4) = -(c1*u(4) + u(3)*(1-u(3)) + a2*u(3)*u(1));
end

%% Boundary Conditions
function res = boundary_conditions(ya, yb, el, er)
% Boundary conditions: connect left equilibrium to right equilibrium
% ya: values at left boundary
% yb: values at right boundary

res = zeros(4, 1);
% At left boundary: approach el (unstable equilibrium)
res(1) = ya(1) - el(1);  % u1(-inf) = 0
res(2) = ya(3) - el(3);  % u2(-inf) = 0

% At right boundary: approach er (stable equilibrium)
res(3) = yb(1) - er(1);  % u1(+inf) = (r+a1)/(r-a1*a2)
res(4) = yb(3) - er(3);  % u2(+inf) = r*(1+a2)/(r-a1*a2)
end
