%% PROJETAR UM CONTROLADOR POR REALIMENTAÇÃO COMPLETA DE ESTADOS DO TIPO PDC INCLUINDO UMA TAXA DE DECAIMENTO

clc; clear;

n = 4; m = 1;
eta = 1e-5;  % taxa de decaimento

% Matrizes A e B (8 regras)
A{1} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 44.7476 13.408 0];    B{1} = [0;1;0;8.9495];
A{2} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 47.2451 9.7444 0];    B{2} = [0;1;0;9.4490];
A{3} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 48.9343 5.9142 0];    B{3} = [0;1;0;9.7869];
A{4} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 49.7861 1.9827 0];    B{4} = [0;1;0;9.9572];
A{5} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 49.7861 -1.9827 0];   B{5} = [0;1;0;9.9572];
A{6} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 48.9343 -5.9142 0];   B{6} = [0;1;0;9.7869];
A{7} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 47.2451 -9.7444 0];   B{7} = [0;1;0;9.4490];
A{8} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 44.7476 -13.408 0];   B{8} = [0;1;0;8.9495];

% Variáveis de decisão
for i = 1:8
    P{i} = sdpvar(n,n,'symmetric');
    Y{i} = sdpvar(m,n,'full');
end

% Restrições LMI
F = [];
for i = 1:8
    F = [F, P{i} >= 1e-6*eye(n)];
    LMI = (A{i}*P{i} + B{i}*Y{i}) + (A{i}*P{i} + B{i}*Y{i})' + 2*eta*P{i};
    F = [F, LMI <= -1e-8*eye(n)];
end

% Objetivo
obj = 0;
for i = 1:8
    obj = obj + trace(P{i});
end

% Solver
options = sdpsettings('solver','sedumi','verbose',1);
sol = optimize(F, obj, options);

if sol.problem == 0
    disp('Sucesso: Controlador TS-PDC com P_i por regra encontrado!');

    for i = 1:8
        K{i} = value(Y{i}) / value(P{i});
        fprintf('\nGanho K da regra %d:\n', i);
        disp(K{i});
    end

    % Simulação
    tspan = [0 10];
    x0 = [0.1; 0; 0.1; 0];
    [t, x] = ode45(@(t,x) pdc_dynamics(x,A,B,K), tspan, x0);

    % Sinal de controle
    u_t = zeros(length(t),1);
    for k = 1:length(t)
        mu_k = fuzzy_membership(x(k,3));
        u_k = 0;
        for j = 1:8
            u_k = u_k + mu_k(j)*(K{j} * x(k,:)');
        end
        u_t(k) = u_k;
    end

    % Gráficos
    figure;
    subplot(5,1,1); plot(t, x(:,1), 'LineWidth', 1.5); ylabel('x_1');
    title('Resposta do sistema com controlador TS-PDC');
    subplot(5,1,2); plot(t, x(:,2), 'LineWidth', 1.5); ylabel('x_2');
    subplot(5,1,3); plot(t, x(:,3), 'LineWidth', 1.5); ylabel('x_3 (ângulo)');
    subplot(5,1,4); plot(t, x(:,4), 'LineWidth', 1.5); ylabel('x_4');
    subplot(5,1,5); plot(t, u_t, 'r', 'LineWidth', 1.5); ylabel('u(t)'); xlabel('Tempo (s)');
    sgtitle('Simulação do sistema TS-PDC com 8 regras');

else
    disp('Problema na otimização:');
    disp(yalmiperror(sol.problem));
    disp(check(F));
end

%% Funções auxiliares 
function mu = fuzzy_membership(x3)
    edges = linspace(-pi, pi, 9);
    mu = zeros(1,8);
    for j = 1:8
        center = (edges(j) + edges(j+1))/2;
        width = edges(j+1) - edges(j);
        mu(j) = max(0, 1 - abs(x3 - center)/(width/2));
    end
    mu = mu / sum(mu + 1e-10); % para normalização
end

function dx = pdc_dynamics(x, A, B, K)
    mu = fuzzy_membership(x(3));
    u = 0;
    for j = 1:8
        u = u + mu(j)*(K{j}*x);
    end
    dx = zeros(length(x),1);
    for j = 1:8
        dx = dx + mu(j)*(A{j}*x + B{j}*u);
    end
end
