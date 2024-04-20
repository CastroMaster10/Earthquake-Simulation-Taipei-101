% Definir parámetros del sistema
M = [6930693 0 0; 0 6930693 0; 0 0 6930693]; % Matriz de masa
K = [-4e11 2e11 0; 2e11 -4e11 2e11; 0 2e11 -2e11]; % Matriz de rigidez
A = inv(M) * K; % Calcular la matriz A

% Calcular los valores propios de la matriz A
[eigenVectors, eigenValues] = eig(A);
% Los valores propios están en la diagonal de eigenValues
omega = sqrt(-diag(eigenValues)); % Frecuencia angular para cada piso
    
%omega = 0.15;
E = 0.5; % Amplitud de la oscilación sísmica
b = [1; 1 ; 1]; % Vector de unos

% Condiciones iniciales (posición y velocidad iniciales)
x0 = [0; 0; 0]; % Posición inicial
v0 = [0; 0; 0]; % Velocidad inicial
initial_conditions = [x0; v0];

% Tiempo de simulación
tspan = [0, 2]; % Simular de 0 a 100 segundos

% Función para las ecuaciones diferenciales
system_odes = @(t, y) [
    y(4:6); % Derivadas de las posiciones (velocidades)
    A*y(1:3) + E * (1) .* cos(1) .* b; % Aceleraciones       
];

% Resolver el sistema de ecuaciones diferenciales
[t, Y] = ode45(system_odes, tspan, initial_conditions);

% Graficar resultados
figure;
plot(t, Y(:, 1), 'r', t, Y(:, 2), 'g', t, Y(:, 3), 'b'); % Posición en función del tiempo
xlabel('Tiempo (s)');
ylabel('Desplazamiento (m)');
title('Respuesta del Sistema SIN Oscilación Sísmica');
legend('Posición Piso 1', 'Posición Piso 2', 'Posición Piso 3');
    
disp(eigenValues)