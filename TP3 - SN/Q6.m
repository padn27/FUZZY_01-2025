%%  OBSERVADOR DE ESTADOS CONSIDERANDO QUE AS VARIÁVEIS PREMISSAS SÃO CONHECIDAS.

clear; clc;

% Parâmetros do sistema
g = 9.81; p = 0.5; c = 1.0; v = 5.0;
theta_max = pi/6;
n = 4;  % número de estados

% Matriz de saída
C = [1 0 0 0];  % Saída medida
epsilon = 1e-3;  % Relaxamento para LMI

% Funções para modelar cos(x3) em 8 setores
cos_vals = cos(linspace(-theta_max, theta_max, 9));
A = cell(8,1);
for i = 1:8
    cos_approx = (cos_vals(i) + cos_vals(i+1)) / 2;
    A{i} = [0 v 0 0;
            0 0 1 0;
            0 0 0 1;
            g/p*cos_approx 0 0 0];
end

% Verificar observabilidade
fprintf('Verificando observabilidade de (A_i, C):\n');
for i = 1:8
    Ob = obsv(A{i}, C);
    r = rank(Ob);
    fprintf('Regra %d: rank = %d\n', i, r);
end

% Variáveis de decisão
P = sdpvar(n,n);      % Matriz de Lyapunov
Y = cell(8,1);         % Ganhos a otimizar
for i = 1:8
    Y{i} = sdpvar(n,1);
end

% Restrição de positividade
Constraints = [P >= epsilon*eye(n)];

% LMIs para cada regra
for i = 1:8
    M = A{i}*P - Y{i}*C;
    Constraints = [Constraints, ...
        M + M' <= -epsilon*eye(n)];
end

% Critério de otimização
Objective = trace(P);  % qualquer função convexa (como norm(Y), etc e etc... e etc.)

% Configurações
ops = sdpsettings('solver','sedumi','verbose',1);
sol = optimize(Constraints, Objective, ops);

% Verificação do resultado
if sol.problem == 0
    disp('Solucionado com sucesso!');
    P_val = value(P);
    L = cell(8,1);
    for i = 1:8
        L{i} = P_val \ value(Y{i});  % L = P⁻¹Y
        fprintf('L{%d} =\n', i);
        disp(L{i});
    end
    
    % Plot dos ganhos L_i
    figure;
    estados = 1:n;
    cores = lines(8); % 8 cores para 8 regras (Bem colorido!)
    hold on;
    for i = 1:8
        plot(estados, L{i}, '-o', 'Color', cores(i,:), 'DisplayName', sprintf('L{%d}', i));
    end
    xlabel('Índice do estado');
    ylabel('Ganho do observador L_i');
    title('Ganhos do Observador para cada regra');
    legend('Location', 'best');
    grid on;
    hold off;
else
    disp('Problema não resolvido:');
    disp(sol.info);
end
