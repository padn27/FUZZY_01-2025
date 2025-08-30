%% CONTROLADOR POR REALIMENTAÇÃO COMPLETA DE ESTADOS DO TIPO U = K.X

clc; clear; close all;

% Parâmetros físicos
p = 0.5; c = 1.0; g = 9.81; v = 5.0;

% Modelo TS: Matrizes A e B (8 regras) 
A{1} = [0 5 0 0;
        0 0 0 0;
        0 0 0 1;
        0 44.7476 13.408 0];
B{1} = [0;1;0;8.9495];

A{2} = [0 5 0 0;
        0 0 0 0;
        0 0 0 1;
        0 47.2451 9.7444 0];
B{2} = [0;1;0;9.4490];

A{3} = [0 5 0 0;
        0 0 0 0;
        0 0 0 1;
        0 48.9343 5.9142 0];
B{3} = [0;1;0;9.7869];

A{4} = [0 5 0 0;
        0 0 0 0;
        0 0 0 1;
        0 49.7861 1.9827 0];
B{4} = [0;1;0;9.9572];

A{5} = [0 5 0 0;
        0 0 0 0;
        0 0 0 1;
        0 49.7861 -1.9827 0];
B{5} = [0;1;0;9.9572];

A{6} = [0 5 0 0;
        0 0 0 0;
        0 0 0 1;
        0 48.9343 -5.9142 0];
B{6} = [0;1;0;9.7869];

A{7} = [0 5 0 0;
        0 0 0 0;
        0 0 0 1;
        0 47.2451 -9.7444 0];
B{7} = [0;1;0;9.4490];

A{8} = [0 5 0 0;
        0 0 0 0;
        0 0 0 1;
        0 44.7476 -13.408 0];
B{8} = [0;1;0;8.9495];

% Ganhos K para cada regra 
K(:,:,1) = [-1.7900 3.2100 4.4034 1.2056];
K(:,:,2) = [-2.4629 2.5371 3.7829 1.2131];
K(:,:,3) = [-4.0581 0.9419 3.2609 1.3342];
K(:,:,4) = [-12.1047 -7.1047 2.8103 2.1195];
K(:,:,5) = [12.1047 17.1047 2.4120 -0.3118];
K(:,:,6) = [4.0581 9.0581 2.0523 0.5050];
K(:,:,7) = [2.4629 7.4629 1.7203 0.6918];
K(:,:,8) = [1.7900 6.7900 1.4070 0.8056];

% Definição dos vértices para cálculo dos pesos (x2 e x3)
x2_vec = [44.7476 47.2451 48.9343 49.7861 49.7861 48.9343 47.2451 44.7476];
x3_vec = [13.4080 9.7444 5.9142 1.9827 -1.9827 -5.9142 -9.7444 -13.4080]; % graus

% Parâmetros da simulação
dt = 0.01;       % passo de tempo (s)
T = 5;           % tempo total (s)
N = T/dt;        % número de passos
x = zeros(4,N);
u = zeros(1,N);

% Condição inicial 
x(:,1) = [0.1; 45; 0; 0]; % [x1, x2, x3 rad, x4]

% Função para calcular pesos fuzzy (coloquei a função auxiliar no final do código)

for k=1:N-1
    % Converter x3 para graus para calcular pesos
    x3_deg = x(3,k)*180/pi;
    
    % Calcula pesos fuzzy
    w = fuzzyWeights(x(2,k), x3_deg, x2_vec, x3_vec);
    
    % Controlador TS: u = - sum(w_i * K_i * x)
    u(k) = 0;
    for i=1:8
        u(k) = u(k) - w(i)*K(:,:,i)*x(:,k);
    end
    
    % Dinâmica ponderada do sistema: A_m e B_m
    A_m = zeros(4);
    B_m = zeros(4,1);
    for i=1:8
        A_m = A_m + w(i)*A{i};
        B_m = B_m + w(i)*B{i};
    end
    
    % Atualiza estado (Euler)
    x(:,k+1) = x(:,k) + dt*(A_m*x(:,k) + B_m*u(k));
end

% Último valor de controle (replicado do penúltimo)
u(N) = u(N-1);

% Plot dos resultados
figure;

subplot(5,1,1);
plot(0:dt:T-dt, x(1,:));
ylabel('x_1 (β) [rad]');
title('Estados do Sistema');

subplot(5,1,2);
plot(0:dt:T-dt, x(2,:));
ylabel('x_2 [rad/s]');

subplot(5,1,3);
plot(0:dt:T-dt, x(3,:));
ylabel('x_3 (α) [rad]');

subplot(5,1,4);
plot(0:dt:T-dt, x(4,:));
ylabel('x_4 [rad/s]');

subplot(5,1,5);
plot(0:dt:T-dt, u);
ylabel('Controle u');
xlabel('Tempo [s]');

sgtitle('Simulação do sistema controlado por realimentação TS');

% --- Função auxiliar para calcular os pesos fuzzy
function w = fuzzyWeights(x2, x3_deg, x2_vec, x3_vec)
    dist = sqrt((x2 - x2_vec).^2 + (x3_deg - x3_vec).^2);
    dist(dist==0) = 1e-6; % para evitar divisão por zero
    w_raw = 1 ./ dist;
    w = w_raw / sum(w_raw); % para normalizar para soma 1
end
