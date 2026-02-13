clc; clear; close all;

% PARAMETERS
L = 1;
alpha_list = [0.01 0.001 0.0001];
T_in = 350; 
T_l = 300; 
T_r = 400;

Nx = 21;
x = linspace(0,L,Nx);
dx = x(2) - x(1);

t_req = [1 5 10 50 100];

% ================= FINITE DIFFERENCE METHODS =================
for a = 1:length(alpha_list)

    alpha = alpha_list(a);

    dt = 0.4*dx^2/alpha;       % stability condition
    Nt = ceil(max(t_req)/dt);
    t = (0:Nt)*dt;
    r = alpha*dt/dx^2;

    % Initial condition
    T_exp = T_in*ones(Nx,1); 
    T_exp([1 end]) = [T_l T_r];

    T_imp = T_exp;
    T_cn  = T_exp;

    % Matrices
    A_imp = diag((1+2*r)*ones(Nx-2,1)) ...
          - diag(r*ones(Nx-3,1),1) ...
          - diag(r*ones(Nx-3,1),-1);

    A_cn = diag((1+r)*ones(Nx-2,1)) ...
         - diag(r/2*ones(Nx-3,1),1) ...
         - diag(r/2*ones(Nx-3,1),-1);

    B_cn = diag((1-r)*ones(Nx-2,1)) ...
         + diag(r/2*ones(Nx-3,1),1) ...
         + diag(r/2*ones(Nx-3,1),-1);

    % Storage
    T_exp_store = zeros(Nx,length(t_req));
    T_imp_store = zeros(Nx,length(t_req));
    T_cn_store  = zeros(Nx,length(t_req));

    t_index = 1;

    for n = 1:length(t)

        % -------- Explicit --------
        Tn = T_exp;
        for i = 2:Nx-1
            Tn(i) = T_exp(i) + r*(T_exp(i+1)-2*T_exp(i)+T_exp(i-1));
        end
        T_exp = Tn;
        T_exp([1 end]) = [T_l T_r];

        % -------- Implicit --------
        b = T_imp(2:end-1);
        b([1 end]) = b([1 end]) + r*[T_l T_r]';
        T_imp(2:end-1) = A_imp\b;

        % -------- Crank-Nicolson --------
        b = B_cn*T_cn(2:end-1);
        b([1 end]) = b([1 end]) + r*[T_l T_r]';
        T_cn(2:end-1) = A_cn\b;

        % -------- Store at requested times --------
        if t_index <= length(t_req) && t(n) >= t_req(t_index)

            T_exp_store(:,t_index) = T_exp;
            T_imp_store(:,t_index) = T_imp;
            T_cn_store(:,t_index)  = T_cn;

            t_index = t_index + 1;
        end
    end

    % -------- Plot --------
    figure('Name',['FDM  alpha = ',num2str(alpha)]);

    for k = 1:length(t_req)
        subplot(1,length(t_req),k); hold on; grid on;
        plot(x,T_exp_store(:,k),'r--','LineWidth',1.5)
        plot(x,T_imp_store(:,k),'b-.','LineWidth',1.5)
        plot(x,T_cn_store(:,k),'k','LineWidth',2)
        title(['t = ',num2str(t_req(k)),' s'])
        xlabel('x'); ylabel('T')
    end
    legend('Explicit','Implicit','Crank-Nicolson')
    sgtitle(['Finite Difference Methods  (\alpha = ',num2str(alpha),')'])
end


% ================= PDEPE SOLUTION =================
nx = 50;
x = linspace(0,L,nx);

for a = 1:length(alpha_list)

    alpha = alpha_list(a);

    % PDE definition
    pdefun = @(x,t,T,dTdx) deal(1, alpha*dTdx, 0);
    icfun  = @(x) T_in;
    bcfun  = @(xl,Tl,xr,Tr,t) deal(Tl-T_l,0,Tr-T_r,0);

    tspan = [0 t_req];
    sol = pdepe(0,pdefun,icfun,bcfun,x,tspan);

    % Plot
    figure('Name',['pdepe  alpha = ',num2str(alpha)]);
    hold on; grid on;

    for k = 1:length(t_req)
        plot(x, sol(k+1,:,1),'LineWidth',2);
    end

    xlabel('x'); ylabel('T')
    title(['pdepe Solution  (\alpha = ',num2str(alpha),')'])
    legend(arrayfun(@(t)['t = ',num2str(t),' s'], t_req,'UniformOutput',false))
end