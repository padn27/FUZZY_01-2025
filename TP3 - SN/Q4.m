%%  SIMULAÇÃO DOS CONTROLADORES CONSIDERANDO UMA TRAJETÓRIA DE REFERÊNCIA E MANTENDO O VEÍCULO EQUILIBRADO NA VERTICAL.

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

% Ganhos K_i encontrados no projeto 
K{1} = 1e4 * [1.2001, -0.4874, -0.7639, -2.9688];
K{2} = 1e4 * [1.1463, -0.7345, -0.6303, -3.4391];
K{3} = 1e4 * [0.9489, -0.8300, -0.4742, -3.9780];
K{4} = 1e4 * [0.4802, -0.7125, -0.3427, -4.4988];
K{5} = 1e4 * [-0.5209, -0.2059, -0.3094, -4.5450];
K{6} = 1e4 * [-1.0333, -0.0399, -0.3861, -4.0629];
K{7} = 1e4 * [-1.2718, -0.1155, -0.4973, -3.5415];
K{8} = 1e4 * [-1.3638, -0.3827, -0.5962, -3.0749];

%% Trajetória de referência para x3 (ângulo)
tspan = [0 10];
ref_fun = @(t) 0.1 * sin(2*pi*0.2*t); % amplitude 0.1 rad, frequência 0.2 Hz

% Condição inicial
x0 = [0; 0; 0; 0]; % Começa equilibrado

%% Simulação
odefun = @(t,x) ts_pdc_ref(t,x,A,B,K,ref_fun);

[t,x] = ode45(odefun, tspan, x0);

% Calcula controle u(t)
u_t = zeros(length(t),1);
ref_t = zeros(length(t),1);
for k=1:length(t)
    mu_k = pertinencia(x(k,3));
    error = x(k,3) - ref_fun(t(k));
    x_err = x(k,:)';
    x_err(3) = error; % substitui x3 pelo erro de referência para controle
    u_val = 0;
    for i=1:length(mu_k)
        u_val = u_val + mu_k(i)*(K{i}*x_err);
    end
    u_t(k) = u_val;
    ref_t(k) = ref_fun(t(k));
end

%% Plot resultados
figure;
subplot(5,1,1);
plot(t,x(:,1)); ylabel('x_1');
title('Estados do sistema com controlador PDC e referência');

subplot(5,1,2);
plot(t,x(:,2)); ylabel('x_2');

subplot(5,1,3);
plot(t,x(:,3)); hold on;
plot(t,ref_t,'r--','LineWidth',1.5);
ylabel('x_3 (ângulo)');
legend('Estado','Referência');

subplot(5,1,4);
plot(t,x(:,4)); ylabel('x_4');

subplot(5,1,5);
plot(t,u_t); ylabel('Controle u'); xlabel('Tempo (s)');

sgtitle('Simulação com referência senoidal e controle PDC');

%% Função pertinência (pesos fuzzy)
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

%% Dinâmica TS com controle PDC e referência
function dx = ts_pdc_ref(t,x,A,B,K,ref_fun)
    mu = pertinencia(x(3));
    error = x(3) - ref_fun(t);
    x_err = x;
    x_err(3) = error; % substitui x3 por erro
    u = 0;
    for i=1:length(mu)
        u = u + mu(i)*(K{i}*x_err);
    end
    dx = zeros(4,1);
    for i=1:length(mu)
        dx = dx + mu(i)*(A{i}*x + B{i}*u);
    end
end
