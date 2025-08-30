%% ESTIMADORES DE ESTADOS ATUANDO SOBRE O SISTEMA CONTROLADO

clear; clc; close all;

% Parâmetros do sistema
g = 9.81; p = 0.5; c = 1.0; v = 5.0;
theta_max = pi/6;
n = 4;    % estados
p_saida = 2; % x1 e x3 medidos
dt = 0.01; % passo de simulação
Tfinal = 5; % segundos

% Matriz de saída
C = [1 0 0 0;
     0 0 1 0];

% Matriz A para as 8 regras
cos_vals = cos(linspace(-theta_max, theta_max, 9));
A = cell(8,1);
for i = 1:8
    cos_approx = (cos_vals(i) + cos_vals(i+1)) / 2;
    A{i} = [0 v 0 0;
            0 0 1 0;
            0 0 0 1;
            g/p*cos_approx 0 0 0];
end

% Matriz B (fixei valor médio ou aproximado de 𝜁42≈9)
B = [0; 1; 0; 9];

% Controle simples (exemplo: u = -Kx)
K = [100, 50, 20, 10];

% Ganhos observador L_i (premissas conhecidas)
L_k = cell(8,1);
L_k{1} = [82.8394    0.0277;
          4.9929   22.2832;
         -0.0214   63.8497;
         19.8525    1.0008];
L_k{2} = [82.8665    0.0277;
          4.9929   22.2832;
         -0.0214   63.8497;
         20.8333    1.0008];
L_k{3} = [82.8849    0.0277;
          4.9929   22.2832;
         -0.0214   63.8497;
         21.4966    1.0008];
L_k{4} = [82.8942    0.0277;
          4.9929   22.2832;
         -0.0214   63.8497;
         21.8312    1.0008];
L_k{5} = [82.8942    0.0277;
          4.9929   22.2832;
         -0.0214   63.8497;
         21.8312    1.0008];
L_k{6} = [82.8849    0.0277;
          4.9929   22.2832;
         -0.0214   63.8497;
         21.4966    1.0008];
L_k{7} = [82.8665    0.0277;
          4.9929   22.2832;
         -0.0214   63.8497;
         20.8333    1.0008];
L_k{8} = [82.8394    0.0277;
          4.9929   22.2832;
         -0.0214   63.8497;
         19.8525    1.0008];

% Ganhos fixos para premissas desconhecidas (média dos L_k)
L_fixed = zeros(n,p_saida);
for i = 1:8
    L_fixed = L_fixed + L_k{i};
end
L_fixed = L_fixed / 8;

% Estado inicial real e estimado
x = zeros(n,1);
x_est_k = zeros(n,1);
x_est_fixed = zeros(n,1);

% Inicializar histórico para plot
time = 0:dt:Tfinal;
X = zeros(n,length(time));
X_est_k = zeros(n,length(time));
X_est_fixed = zeros(n,length(time));

for k = 1:length(time)
    % Premissa real: ângulo theta = x(3)
    theta = x(3);
    
    % Índice da regra ativa (1 a 8)
    sector_width = 2*theta_max/8;
    i = min(max(floor((theta + theta_max)/sector_width) + 1,1),8);
    
    % Controle
    u = -K * x;
    
    % Dinâmica do sistema
    dx = A{i}*x + B*u;
    
    % Atualiza estado real (Euler)
    x = x + dt*dx;
    
    % Saída com ruído
    y = C * x + 0.01*randn(p_saida,1);
    
    % Estimador premissa conhecida
    dx_est_k = A{i}*x_est_k + B*u + L_k{i}*(y - C*x_est_k);
    x_est_k = x_est_k + dt*dx_est_k;
    
    % Estimador premissa desconhecida (ganho fixo)
    dx_est_fixed = A{i}*x_est_fixed + B*u + L_fixed*(y - C*x_est_fixed);
    x_est_fixed = x_est_fixed + dt*dx_est_fixed;
    
    % Armazena para plot
    X(:,k) = x;
    X_est_k(:,k) = x_est_k;
    X_est_fixed(:,k) = x_est_fixed;
end

% Plot resultados
figure;
for idx = 1:n
    subplot(n,1,idx);
    hold on;
    plot(time, X(idx,:), 'k', 'LineWidth', 1.5);
    plot(time, X_est_k(idx,:), 'b--', 'LineWidth', 1);
    plot(time, X_est_fixed(idx,:), 'r:', 'LineWidth', 1);
    ylabel(sprintf('x_%d', idx));
    if idx == 1
        legend('Real','Estimador Premissa Conhecida','Estimador Premissa Desconhecida');
        title('Estados e Estimativas');
    end
    grid on;
end
xlabel('Tempo (s)');
