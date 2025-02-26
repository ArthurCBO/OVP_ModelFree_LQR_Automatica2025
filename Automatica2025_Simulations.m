clc
clear all
close all

%%
% General System Setup

n = 5;
m = 3;
k = 10;

Q = 1e0*eye(n);
R = 1e0*eye(m);

% A = randHurwitz(n);
A = [-5.2373   -0.3452   -0.6653   -0.6715   -0.3288
     -0.3452   -5.4889   -0.8060   -0.3889   -0.5584
     -0.6653   -0.8060   -5.0377   -0.5735   -0.5100
     -0.6715   -0.3889   -0.5735   -5.3354   -0.6667
     -0.3288   -0.5584   -0.5100   -0.6667   -5.4942 ];
B = [zeros(n-m, m) ; eye(m)];
C = eye(n);
D = [];

[Kstar, ~, ~] = lqr(A,B,Q,R);
Kstar = -Kstar;

[PsiStar, SigStar, PhiStar] = svd(Kstar);
GammaBal = RandOrthMat(k);

probParams = [];
probParams.A = A;
probParams.B = B;
probParams.R = R;
probParams.Q = Q;
probParams.m = m;
probParams.n = n;
probParams.k = k;

%%
% General parameters definition

islog = true;

if islog
    RangeofImb = [0,2,100];
    l = logspace(RangeofImb(1), RangeofImb(2), RangeofImb(3));
else
    RangeofImb = [1,100,100];
    l = linspace(RangeofImb(1), RangeofImb(2), RangeofImb(3));
end
rangel = length(l);

%%
% Simulations for eta>1

eta = 5;

if islog
    tsim = logspace(-6, 2, 1000);
else
    tsim = linspace(0,100,1000);
end

SigK = sqrt(SigStar*abs(eta));
K20_balanced = PsiStar*[SigK, zeros(m, k-n)]*GammaBal.';
K10_balanced = sign(eta)*GammaBal*[SigK ; zeros(k-m, n)]*PhiStar.';


[K10,K20] = ComputeInitializations(K10_balanced,K20_balanced,RangeofImb,islog);

[K,J] = SimulateRegularCase(eta*Kstar,tsim,probParams);
[K1,K2,Jovp] = SimulateOvpCase(K10,K20,tsim,probParams);

ttl = ['\eta=' num2str(eta)];
plotFigures('eta5', tsim, J, Jovp, l, ttl)

%%%%%%%%%%%
% 

eta = 20;

if islog
    tsim = logspace(-6, 2, 1000);
else
    tsim = linspace(0,100,1000);
end

SigK = sqrt(SigStar*abs(eta));
K20_balanced = PsiStar*[SigK, zeros(m, k-n)]*GammaBal.';
K10_balanced = sign(eta)*GammaBal*[SigK ; zeros(k-m, n)]*PhiStar.';


[K10,K20] = ComputeInitializations(K10_balanced,K20_balanced,RangeofImb,islog);

[K,J] = SimulateRegularCase(eta*Kstar,tsim,probParams);
[K1,K2,Jovp] = SimulateOvpCase(K10,K20,tsim,probParams);

ttl = ['\eta=' num2str(eta)];
plotFigures('eta20', tsim, J, Jovp, l, ttl)

%%
% Simulations for 1>eta>0

eta = 0.9;

if islog
    tsim = logspace(log10(5e-5), log10(500), 1000);
else
    tsim = linspace(0,100,1000);
end

SigK = sqrt(SigStar*abs(eta));
K20_balanced = PsiStar*[SigK, zeros(m, k-n)]*GammaBal.';
K10_balanced = sign(eta)*GammaBal*[SigK ; zeros(k-m, n)]*PhiStar.';


[K10,K20] = ComputeInitializations(K10_balanced,K20_balanced,RangeofImb,islog);

[K,J] = SimulateRegularCase(eta*Kstar,tsim,probParams);
[K1,K2,Jovp] = SimulateOvpCase(K10,K20,tsim,probParams);

ttl = ['\eta=' num2str(eta)];
plotFigures('eta09', tsim, J, Jovp, l, ttl)

%%%%%%%%%%%%%%%%%%
% 

eta = 0.1;

if islog
    tsim = logspace(-4, 3, 1000);
else
    tsim = linspace(0,100,1000);
end

SigK = sqrt(SigStar*abs(eta));
K20_balanced = PsiStar*[SigK, zeros(m, k-n)]*GammaBal.';
K10_balanced = sign(eta)*GammaBal*[SigK ; zeros(k-m, n)]*PhiStar.';


[K10,K20] = ComputeInitializations(K10_balanced,K20_balanced,RangeofImb,islog);

[K,J] = SimulateRegularCase(eta*Kstar,tsim,probParams);
[K1,K2,Jovp] = SimulateOvpCase(K10,K20,tsim,probParams);

ttl = ['\eta=' num2str(eta)];
plotFigures('eta01', tsim, J, Jovp, l, ttl)

%%
% Simulations for 1>eta>0

eta = -0.1;

if islog
    tsim = logspace(-4, 4, 1000);
else
    tsim = linspace(0,1000,1000);
end

SigK = sqrt(SigStar*abs(eta));
K20_balanced = PsiStar*[SigK, zeros(m, k-n)]*GammaBal.';
K10_balanced = sign(eta)*GammaBal*[SigK ; zeros(k-m, n)]*PhiStar.';


[K10,K20] = ComputeInitializations(K10_balanced,K20_balanced,RangeofImb,islog);

[K,J] = SimulateRegularCase(eta*Kstar,tsim,probParams);
[K1,K2,Jovp] = SimulateOvpCase(K10,K20,tsim,probParams);

ttl = ['\eta=' num2str(eta)];
plotFigures('etam01', tsim, J, Jovp, l, ttl)

%%
% 

eta = -20;

if islog
    tsim = logspace(-7, 3, 1000);
else
    tsim = linspace(0,100,1000);
end

SigK = sqrt(SigStar*abs(eta));
K20_balanced = PsiStar*[SigK, zeros(m, k-n)]*GammaBal.';
K10_balanced = sign(eta)*GammaBal*[SigK ; zeros(k-m, n)]*PhiStar.';


[K10,K20] = ComputeInitializations(K10_balanced,K20_balanced,RangeofImb,islog);

[K,J] = SimulateRegularCase(eta*Kstar,tsim,probParams);
[K1,K2,Jovp] = SimulateOvpCase(K10,K20,tsim,probParams);
% 
ttl = ['\eta=' num2str(eta)];
plotFigures('etam20', tsim, J, Jovp, l, ttl)

%%
% Investigate why the solution does not approach 0

c=arrayfun(@(i)trace(K1(:,:,i,1)*K1(:,:,i,1)'-K2(:,:,i,1)'*K2(:,:,i,1)), 1:1000);
figure('Renderer', 'painters')
f = gcf;
plot(tsim, c, 'LineWidth',2);
xlabel('Training Time');
ylabel('Trace of Invariant (Tr(C))');
title(ttl);
set(gca, 'XScale', 'log', 'YScale', 'linear');
figname = [date '_InvariantDuringTraining_' 'etam20' '.png'];
exportgraphics(f,figname,'Resolution',600)

%%
% 

k_ = pagenorm(pagemtimes(K2(:,:,:,1),K1(:,:,:,1)), 'fro');
ks = norm(Kstar,'fro');
figure('Renderer', 'painters', 'Position', [680, 558, 420, 315])
f = gcf;
hold on
plot(tsim, squeeze(k_), 'LineWidth',2);
plot(tsim, ks*ones(length(tsim),1), 'k--', 'LineWidth', 2);
xlabel('Training Time');
ylabel('||K_2(t)K_1(t)||_F');
title(ttl);
set(gca, 'XScale', 'log', 'YScale', 'linear');
box on
figname = [date '_NormProdDuringTraining_' 'etam20' '.png'];
exportgraphics(f,figname,'Resolution',600)

disp(Jfunc(A, B, K2(:,:,end,1)*K1(:,:,end,1),R, Q))
disp(Jfunc(A, B, Kstar, R, Q))
disp(Jfunc(A, B, zeros(size(Kstar)), R, Q))

%%
% Simulations for Different Initialization

eta = diag([20 ; 0.1 ; -20]);

if islog
    tsim = logspace(-6, log10(500), 1000);
else
    tsim = linspace(0,1000,1000);
end

SigK = sqrt(abs(eta)*SigStar);
K20_balanced = PsiStar*[SigK, zeros(m, k-n)]*GammaBal.';
K10_balanced = GammaBal*[sign(eta)*SigK ; zeros(k-m, n)]*PhiStar.';


[K10,K20] = ComputeInitializations(K10_balanced,K20_balanced,RangeofImb,islog);

[K,J] = SimulateRegularCase(K20_balanced*K10_balanced,tsim,probParams);
[K1,K2,Jovp] = SimulateOvpCase(K10,K20,tsim,probParams);

ttl = '\eta=diag([20, 0.1, -20])';
plotFigures('eta2001m20', tsim, J, Jovp, l, ttl)

%%
% General Functions Definitions

function [K10,K20] = ComputeInitializations(K10_balanced,K20_balanced,RangeofImb,islog)

    K20 = [];
    K10 = [];

    if islog
        l = logspace(RangeofImb(1), RangeofImb(2), RangeofImb(3));
    else
        l = linspace(RangeofImb(1), RangeofImb(2), RangeofImb(3));
    end

    for ind = 1:RangeofImb(3)
        K20(:,:,ind) = K20_balanced/l(ind);
        K10(:,:,ind) = K10_balanced*l(ind);
    end

end

function [K,J] = SimulateRegularCase(K0, tsim, probParams)

    k0 = K0(:);
    [~, k] = ode45(@(t,x)kdot_func(x, probParams.A, probParams.B, probParams.R, probParams.Q, probParams.m, probParams.n, probParams.k), tsim, k0);
    
    J = [];
    K = [];
    for ind = 1:length(tsim)
        K(:,:,ind) = reshape(k(ind, :), [probParams.m,probParams.n]);
        J(ind) = Jfunc(probParams.A,probParams.B,K(:,:,ind),probParams.R,probParams.Q);
    end

end

function [K1, K2, Jovp] = SimulateOvpCase(K10, K20, tsim, probParams)
    rangel = size(K10, 3);
    K1 = [];
    K2 = [];
    Jovp = [];

    for ind = 1:rangel
        K10_ = K10(:,:,ind);
        K20_ = K20(:,:,ind);
        k0unb = [K10_(:) ; K20_(:)];
        [~, kunb_] = ode15s(@(t,x)kdot_func_ovp(x, probParams.A, probParams.B, probParams.R,probParams.Q, probParams.m, probParams.n, probParams.k), tsim, k0unb);
        kunb(:,:,ind) = kunb_;

        for indj = 1:length(tsim)
        
            K1(:,:,indj,ind) = reshape(kunb(indj, 1:probParams.n*probParams.k, ind), [probParams.k,probParams.n]);
            K2(:,:,indj,ind) = reshape(kunb(indj, probParams.n*probParams.k+1:end, ind), [probParams.m,probParams.k]);
        
            Jovp(indj, ind) = Jfunc(probParams.A,probParams.B,K2(:,:,indj,ind)*K1(:,:,indj,ind),probParams.R,probParams.Q);
        
        end
    end
end

%%
% System's Dynamics

function dkdt = kdot_func_ovp(k_,A,B,R,Q,m,n,k)
    
    k1 = k_(1:n*k);
    k2 = k_(n*k+1:end);
    K1 = reshape(k1,[k,n]);
    K2 = reshape(k2,[m,k]);
    
    Kbar = K2*K1;
    [PK,~,~] = icare(A+B*Kbar, zeros(n,m), Q, -inv(R), Kbar.', eye(n), zeros(n));
    [L,~,~] = icare((A+B*Kbar).', zeros(n,m), eye(n), -eye(m), zeros(n,m), eye(n), zeros(n));
    
    try
        K1dot = -2*((B*K2).'*PK+K2.'*R*K2*K1)*L;
        K2dot = -2*(B.'*PK+R*K2*K1)*L*K1.';
    catch
        disp('Something went wrong during simulation. Probably unstable initialization.');
    end

    dkdt = [K1dot(:) ; K2dot(:)];
    
end

function dkdt = kdot_func(k_,A,B,R,Q,m,n,k)
    K = reshape(k_, [m,n]);

    At = A+B*K;
    Bt = zeros(n,m);
    Qt = Q;
    Rt = -inv(R);
    St = K.';
    Et = eye(n);
    Gt = zeros(n);
    [Pk,~,~] = icare(At,Bt,Qt,Rt,St,Et,Gt);

    At = (A+B*K).';
    Bt = zeros(n,m);
    Qt = eye(n);
    Rt = eye(m);
    St = zeros(n,m);
    Et = eye(n);
    Gt = zeros(n);
    [Lk,~,~] = icare(At,Bt,Qt,Rt,St,Et,Gt);
    
    try
        Kdot = -2*(B.'*Pk+R*K)*Lk;
    catch 
        disp('Something went wrong during simulation. Probably unstable initialization.');
    end
    dkdt = Kdot(:);
end

function J = Jfunc(A,B,K,R,Q)
    sizeB = size(B);
    n = sizeB(1);
    m = sizeB(2);
    
    At = A+B*K;
    Bt = zeros(n,m);
    Qt = Q;
    Rt = -inv(R);
    St = K.';
    Et = eye(n);
    Gt = zeros(n);
    [Pk,~,~] = icare(At,Bt,Qt,Rt,St,Et,Gt);
    J = trace(Pk);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot and save figures

function plotFigures(prefix, t, J, Jovp, l, ttl)

    rangel = length(l);

    %%%
    % Color gradient definition for plots
    
    lightBLUE = [0.356862745098039,0.811764705882353,0.956862745098039];
    darkBLUE = [0.0196078431372549,0.0745098039215686,0.670588235294118];
     
    blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
    
    map = [];
    c=0.25;
    lcut = unique(round(c*l)/c);
    rangelcut = length(lcut);
    
    
    for ind =1:rangelcut
        map = [map ; blueGRADIENTflexible(ind, rangelcut)];
    end
    
    close all
    % figure('Renderer', 'painters')
    % clf
    % hold on
    % % plot(t, distbal, 'LineWidth',2);
    % % leg = [];
    % clear leg
    % for ind = 1:rangel
    %     plot(t, Jovp(:, ind), 'color', blueGRADIENTflexible(ind,rangel), 'LineWidth',2);
    % %     legend(['Overpar, Imbalance = ' num2str(l(ind))]);
    %     leg{ind} = ['Overparametrized, \mu = ' num2str(l(ind))];
    % end
    % 
    % plot(t, J, 'r', 'LineWidth',2);
    % leg{ind+1} = ['Not Overparametrized'];
    % xlabel('Training Time')
    % ylabel('Cost J(K)')
    % title(ttl);
    % colormap(map)
    % cb = colorbar;
    % set(gca, 'clim', [0.5 rangelcut+0.5]);
    % set(cb, 'ticks', 1:rangelcut, 'ticklabels', lcut);
    % colorTitleHandle = get(cb,'Title');
    % titleString = {['  Level of\newlineImbalance']};
    % set(colorTitleHandle ,'String',titleString);
    % box on
    % 
    % f = gcf;
    % figname = [date '_CostvsTraining_' prefix '.png'];
    % exportgraphics(f,figname,'Resolution',600)
    
    figure('Renderer', 'painters', 'Position', [680, 558, 420, 315]);
    % figure(2)
    clf
    hold on
    % plot(t, distbal, 'LineWidth',2);
    % leg = [];
    clear leg2
    for ind = 1:rangel
        plot(t, Jovp(:, ind), 'color', blueGRADIENTflexible(ind,rangel), 'LineWidth',2);
        leg2{ind} = ['Overparametrized, \mu = ' num2str(l(ind))];
    end
    
    plot(t, J(:), 'r', 'LineWidth',2);
    leg2{ind+1} = ['Not Overparametrized'];
    xlabel('Training Time')
    ylabel('Cost J(K)')
    title(ttl);
    
    colormap(map)
    cb = colorbar;
    set(gca, 'clim', [0.5 rangelcut+0.5]);
    set(cb, 'ticks', 1:rangelcut, 'ticklabels', lcut);
    colorTitleHandle = get(cb,'Title');
    titleString = {['  Level of\newlineImbalance']};
    set(colorTitleHandle ,'String',titleString);
    
    set(gca, 'XScale', 'log', 'YScale', 'linear')
    xlim([t(1), t(end)]);
    
    box on
    f = gcf;
    % annotation('textarrow', [0.66, 0.2], [0.66, 0.66], 'LineWidth', 5, 'HeadLength', 20, 'HeadWidth', 20, 'String',{' Increasing', ' Imbalance'});
    figname = [date '_CostvsTrainingLog_' prefix '.png'];
    % drawnow();
    % exportgraphics(f,figname,'Resolution',600)
end