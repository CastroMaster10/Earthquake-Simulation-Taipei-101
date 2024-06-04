clear; close all; clc;

n = 101 ;
k = 8e5;
c = 4e3;
m = 6937.4256;
M = diag(ones(1, n) * m);
K = zeros(n);
K(1:n+1:end) = -k * 2;


for i = 1:n-1
    K(i, i+1) = k; % Asignar a la derecha
    K(i+1, i) = k; % Asignar a la izquierda
end
K(n, n) = -k;
C = zeros(n);
C(1:n+1:end) = -c * 2;
C(1, 1) = -c;
for i = 1:n-1
    C(i, i+1) = c; % Asignar a la derecha
    C(i+1, i) = c; % Asignar a la izquierda
end
C(n, n) = -c;
b = zeros(n,1);
b(1) = 1;

fb = ones(n*2,1);
% State-space matrices
A = [zeros(n), eye(n); inv(M)*K, inv(M)*C];
B = [zeros(n, 1); inv(M) * b];
C = [eye(n),zeros(n)];
D = zeros(n, 1);

% Calculate omega and E
omega = sqrt(k / m);
E = 0.1; % Amplitude of the seismic oscillation

% Initial conditions (including earthquake impact)
initial_displacement = E * cos(omega* 0); % Displacement at t=0
initial_velocity = -E * omega * sin(omega * 0); % Velocity at t=0

x0 = initial_displacement * b; % Initial displacement vector
v0 = initial_velocity * b; % Initial velocity vector

initial_conditions = [x0;v0];

% Save matrices to workspace for Simulink
save('state_space_matrices.mat', 'A', 'B', 'C', 'D', 'initial_conditions');

