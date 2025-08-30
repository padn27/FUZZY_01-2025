%% ESTIMATIVA DOS ESTADOS PELO OBSERVADOR TS.

clear; clc;

% Parâmetros do sistema
g = 9.81; p = 0.5; c = 1.0; v = 5.0;
theta_max = pi/6;
n = 4;  % número de estados
m = 1;  % dimensão da saída
r = 8;  % número de regras TS

% Matrizes A e B (8 regras) 
A = cell(r,1);
B = cell(r,1);

A{1} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 44.7476 13.408 0];    B{1} = [0;1;0;8.9495];
A{2} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 47.2451 9.7444 0];    B{2} = [0;1;0;9.4490];
A{3} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 48.9343 5.9142 0];    B{3} = [0;1;0;9.7869];
A{4} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 49.7861 1.9827 0];    B{4} = [0;1;0;9.9572];
A{5} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 49.7861 -1.9827 0];   B{5} = [0;1;0;9.9572];
A{6} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 48.9343 -5.9142 0];   B{6} = [0;1;0;9.7869];
A{7} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 47.2451 -9.7444 0];   B{7} = [0;1;0;9.4490];
A{8} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 44.7476 -13.408 0];   B{8} = [0;1;0;8.9495];

% Matriz de saída
C = [1 0 0 0];

% Ganhos L_i (obtidos do solver)
L = cell(r,1);
L{1} = [80.5658; 13.0140; 15.5003; 14.5842];
L{2} = [80.6187; 13.3124; 15.7703; 15.1651];
L{3} = [80.6535; 13.5142; 15.9529; 15.5579];
L{4} = [80.6719; 13.6159; 16.0450; 15.7561];
L{5} = [80.6719; 13.6159; 16.0450; 15.7561];
L{6} = [80.6535; 13.5142; 15.9529; 15.5579];
L{7} = [80.6187; 13.3124; 15.7703; 15.1651];
L{8} = [80.5658; 13.0140; 15.5003; 14.5842];

% Função para calcular pertinências mu_i baseado em theta = x(3)
% Dividi o intervalo [-theta_max, theta_max] em 8 setores com função triangular
function mu = calcula_mu(theta, theta_max)
    % Pontos de quebra para 9 limites (8 setores)
    limites = linspace(-theta_max, theta_max, 9);
    mu = zeros(8,1);
    for i = 1:8
        left = limites(i);
        center = (limites(i)+limites(i+1))/2;
        right = limites(i+1);
        if theta >= left && theta <= center
            mu(i) = (theta - left)/(center - left);
        elseif theta > center && theta <= right
            mu(i) = (right - theta)/(right - center);
        else
            mu(i) = 0;
        end
    end
    % Normaliza para somar 1
    mu = mu / sum(mu);
end

% Simulação do observador TS
dt = 0.001;       % passo de integração
T = 5;            % tempo total
N = T/dt;         % número de passos
t = linspace(0, T, N);

% Inicialização do estado estimado
x_hat = zeros(n,1);

% Estado real do sistema 
x_real = [0; 0; 0; 0]; 

% Entrada constante para teste
u = 1;

% Vetores para armazenar dados
X_hat = zeros(n,N);
X_real = zeros(n,N);
Y = zeros(m,N);

for k=1:N
    theta = x_hat(3); % variável premissa do observador
    
    % Calcula pertinências
    mu = calcula_mu(theta, theta_max);
    
    % Mede saída (com ruído zero)
    y = C * x_real;
    
    % Calcula a derivada do estado estimado do observador
    x_hat_dot = zeros(n,1);
    for i=1:r
        x_hat_dot = x_hat_dot + mu(i) * ( A{i}*x_hat + B{i}*u + L{i}*(y - C*x_hat) );
    end
    
    % Atualiza o estado estimado (Euler)
    x_hat = x_hat + dt * x_hat_dot;
    
    % Atualiza o estado real do sistema para simulação (considerando modelo nominal)
    % Para simplicidade A{4} e B{4} como sistema real
    x_real_dot = A{4}*x_real + B{4}*u;
    x_real = x_real + dt * x_real_dot;
    
    % Armazena dados
    X_hat(:,k) = x_hat;
    X_real(:,k) = x_real;
    Y(:,k) = y;
end

% Plota resultados
figure;
for i=1:n
    subplot(n,1,i)
    plot(t,X_real(i,:), 'b', 'LineWidth',1.5); hold on
    plot(t,X_hat(i,:), 'r--', 'LineWidth',1.5);
    ylabel(sprintf('x_%d', i))
    if i==1
        legend('Estado real','Estimado')
    end
end
xlabel('Tempo (s)')
sgtitle('Estimativa dos Estados pelo Observador TS')

