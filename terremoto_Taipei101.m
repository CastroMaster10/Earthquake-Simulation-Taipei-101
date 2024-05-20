clear; close all; clc;

n = 101;

k = 8e5;
c= 4e3;
m = 6937.4256;


M = diag(ones(1, n) * m);
K = zeros(n);
K(1:n+1:end) = -k * 2;

for i = 1:n-1
    K(i, i+1) = k; % Asignar a la derecha
    K(i+1, i) = k; % Asignar a la izquierda
end

K(n, n) = -k;

A = inv(M) * K; % Calcular la matriz A

C = zeros(n);

C(1:n+1:end) = -c * 2;
C(1, 1) = -c;
for i = 1:n-1
    C(i, i+1) = c; % Asignar a la derecha
    C(i+1, i) = c; % Asignar a la izquierda
end

C(n, n) = -c;

D = inv(M) * C;

% Calcular los valores propios de la matriz A
[eigenVectors, eigenValues] = eig(A);
% Los valores propios están en la diagonal de eigenValues
%omega = sqrt(-diag(eigenValues)); % Frecuencia angular para cada piso
    
omega = sqrt(k/m);
E = 0.1; % Amplitud de la oscilación sísmica
b = zeros(n, 1); % Vector de unos 
b(1)=1;

% Condiciones iniciales (posición y velocidad iniciales)
x0 = zeros(n, 1); % Posición inicial
v0 = zeros(n, 1); % Velocidad inicial
initial_conditions = [x0; v0];

% Tiempo de simulación
tspan = [0, 60]; % Simular de 0 a 100 segundos

% Función para las ecuaciones diferenciales
system_odes = @(t, y) [
    y(n+1:end); % Derivadas de las posiciones (velocidades)
    A*y(1:n) + (t > 0 & t < 4).* E * (omega.^2) .* cos(omega * t) .* b + D * y(n+1:end); % Aceleraciones       
];

% Resolver el sistema de ecuaciones diferenciales
[t, Y] = ode45(system_odes, tspan, initial_conditions);

% Crear la leyenda
leyenda = cell(n, 1);
for i = 1:n
    leyenda{i} = ['Posición Piso ' num2str(i)];
end

% Graficar resultados
figure;
hold on;
for i = 1:n
    plot(t, Y(:, i), 'Color', rand(1,3));
end
hold off;
xlabel('Tiempo (s)');
ylabel('Desplazamiento (m)');
title('Respuesta del Sistema a una Oscilación Sísmica');
legend(leyenda);

disp(omega);