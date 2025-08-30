%% SIMULAÇÃO DOS CONTROLADORES CONSIDERANDO UMA TRAJETÓRIA DE REFERÊNCIA E CONSIDERANDO TRAJETÓRIA PARA X3D = K SIN(ꞶT) E X4D = K COS(ꞶT). 
%% QUAIS OS VALORES MÁXIMOS DE K E Ꞷ O SISTEMA CONSEGUE RASTREAR.

clc; clear; close all;

%% Dados do sistema (copiados do projeto do controlador)
n = 4; m = 1;

% Matrizes A e B (8 regras)
A{1} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 44.7476 13.408 0];    B{1} = [0;1;0;8.9495];
A{2} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 47.2451 9.7444 0];    B{2} = [0;1;0;9.4490];
A{3} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 48.9343 5.9142 0];    B{3} = [0;1;0;9.7869];
A{4} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 49.7861 1.9827 0];    B{4} = [0;1;0;9.9572];
A{5} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 49.7861 -1.9827 0];   B{5} = [0;1;0;9.9572];
A{6} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 48.9343 -5.9142 0];   B{6} = [0;1;0;9.7869];
A{7} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 47.2451 -9.7444 0];   B{7} = [0;1;0;9.4490];
A{8} = [0 5 0 0; 0 0 0 0; 0 0 0 1; 0 44.7476 -13.408 0];   B{8} = [0;1;0;8.9495];

% Ganhos K_i 
K{1} = 1e4 * [1.2001, -0.4874, -0.7639, -2.9688];
K{2} = 1e4 * [1.1463, -0.7345, -0.6303, -3.4391];
K{3} = 1e4 * [0.9489, -0.8300, -0.4742, -3.9780];
K{4} = 1e4 * [0.4802, -0.7125, -0.3427, -4.4988];
K{5} = 1e4 * [-0.5209, -0.2059, -0.3094, -4.5450];
K{6} = 1e4 * [-1.0333, -0.0399, -0.3861, -4.0629];
K{7} = 1e4 * [-1.2718, -0.1155, -0.4973, -3.5415];
K{8} = 1e4 * [-1.3638, -0.3827, -0.5962, -3.0749];

%% Referência
K_ref = 0.1;
omega_ref = 1.0;
tspan = [0 15];
x0 = [0; 0; 0; 0];

%% Simulação
[t,x] = ode45(@(t,x) ts_pdc_param_ref(t,x,A,B,K,K_ref,omega_ref), tspan, x0);

% Referência
x3d = K_ref * sin(omega_ref*t);
x4d = K_ref * cos(omega_ref*t);

% Sinal de controle
u_t = zeros(length(t),1);
for k=1:length(t)
    mu_k = pertinencia(x(k,3));
    error = [x(k,3)-x3d(k); x(k,4)-x4d(k)];
    x_err = x(k,:)';
    x_err(3) = error(1);
    x_err(4) = error(2);
    for i=1:length(mu_k)
        u_t(k) = u_t(k) + mu_k(i)*(K{i}*x_err);
    end
end

%% Gráficos principais
figure;
subplot(5,1,1); plot(t,x(:,1)); ylabel('x_1');
title('Estados do sistema com controlador PDC e referência paramétrica');

subplot(5,1,2); plot(t,x(:,2)); ylabel('x_2');

subplot(5,1,3); plot(t,x(:,3), 'b'); hold on;
plot(t,x3d,'r--','LineWidth',1.5); ylabel('x_3 (ângulo)');
legend('Estado','Referência');

subplot(5,1,4); plot(t,x(:,4), 'b'); hold on;
plot(t,x4d,'r--','LineWidth',1.5); ylabel('x_4');
legend('Estado','Referência');

subplot(5,1,5); plot(t,u_t); ylabel('u(t)'); xlabel('Tempo (s)');

sgtitle(sprintf('Simulação PDC com referência: K=%.2f, ω=%.2f', K_ref, omega_ref));

%% Plano de fase x3 vs x4
figure;
plot(x(:,3), x(:,4), 'b-', 'LineWidth', 1.5); hold on;
plot(x3d, x4d, 'r--', 'LineWidth', 1.5);
xlabel('x_3 (ângulo)'); ylabel('x_4 (vel. angular)');
legend('Trajetória real', 'Referência');
title('Plano de fase: x_3 vs x_4'); grid on;

%% Função pertinência
function mu = pertinencia(x3)
    n_rules = 8;
    edges = linspace(-pi, pi, n_rules+1);
    mu = zeros(1,n_rules);
    for i=1:n_rules
        center = (edges(i)+edges(i+1))/2;
        width = (edges(i+1)-edges(i));
        mu(i) = max(0, 1 - abs(x3 - center)/(width/2));
    end
    s = sum(mu);
    if s > 0
        mu = mu/s;
    else
        mu = ones(1,n_rules)/n_rules;
    end
end

%% Dinâmica do sistema
function dx = ts_pdc_param_ref(t,x,A,B,K,K_ref,omega_ref)
    mu = pertinencia(x(3));
    x3d = K_ref * sin(omega_ref*t);
    x4d = K_ref * cos(omega_ref*t);
    error = [x(3)-x3d; x(4)-x4d];
    x_err = x;
    x_err(3) = error(1);
    x_err(4) = error(2);
    u = 0;
    for i=1:length(mu)
        u = u + mu(i)*(K{i}*x_err);
    end
    dx = zeros(4,1);
    for i=1:length(mu)
        dx = dx + mu(i)*(A{i}*x + B{i}*u);
    end
end
