%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2019 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear;
close all;
%
Flag.save = 0;
%
global sqrtK til_f matSigma bar_pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
param.nominal_load = 10.0;
param.alpha = 5.0;
param.vol = 1000; % *(10^5) [mm^3]
param.num_Pareto = 30;
%
eeee = 2.0;
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data of the structure
% --->
nx = 4;
ny = 4;

filename.std = strcat('x',num2str(nx),'_y',num2str(ny));

[dll,matH,coord_x,ir,irr,ird] = member(nx, ny, 1.5);

nk = size(coord_x,1); num.node   = nk;
nd = size(matH,1);    num.degree = nd;
nm = size(matH,2);    num.member = nm;
% <---
% Data of the structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vector
% --->
Idx.deg_of_load = 2 * (nx*(ceil((ny+1)/2)));

til_f = sparse(zeros(nd,1));
til_f(Idx.deg_of_load) = -param.nominal_load;
%
nn = length(Idx.deg_of_load);
matM = sparse( (Idx.deg_of_load-1), (1:nn), ones(1,nn), nd, nn);
matM = [matM, 0.2 * sparse( Idx.deg_of_load, (1:nn), ones(1,nn), nd, nn)];
matSigma = ((param.alpha/3)^2) * matM * (matM');
clear nn
% <---
% Load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness matrix
% --->
sqrtK = zeros(nd,nm);
for i=1:nm
    sqrtK(:,i) = sqrt(eeee / dll(i)) * matH(:,i);
end
sqrtK = sparse(sqrtK);
%
vec_cs = 20 * ones(nm,1);
matK   = sqrtK * diag(vec_cs) * sqrtK';
vec_u  = matK \ til_f;
vec_sigma = (eeee ./ dll) .* (matH' * vec_u);
fprintf(' ============================================= \n');
fprintf('   max |f| = %3.3f [kN]\n',...
    norm(til_f,inf) * 10);
fprintf('   compliance of init.sol. = %3.3f [J]  \n',...
    (til_f' * vec_u * 10) );
% <---
% Stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_aMM = (param.alpha^2) * matM * (matM');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nominal optimization with mean load
% --->
fprintf(' ==  Nominal optimization w/ load at center of ellipsoid >>> \n');
cvx_begin sdp quiet
    cvx_solver sedumi
    cvx_precision best
    variable nom_x(nm,1)
    variable nom_w(1,1)
    minimize( nom_w )
    subject to
        [nom_w, til_f';
            til_f, (sqrtK * diag(nom_x) * sqrtK')] >= 0;
        dll' * nom_x <= param.vol;
        nom_x >= 0;
cvx_end

nom_matK = sqrtK * diag(nom_x) * sqrtK';
[nom_w_worst,~] = worst_analysis(nom_matK, mat_aMM, til_f);

min_pi = til_f' * (nom_matK \ til_f);
fprintf('     compliance for mean force = %8.3f [J] \n',...
    (min_pi * 10) );
fprintf('     worst-case compliance     = %8.3f [J] \n',...
    (nom_w_worst * 10) );
fprintf('     variance of compliance = %3.3f [J^2] \n',...
    obj_variance(nom_x) * 10^2 );
fprintf(' ============================================= \n');
% <---
% Nominal optimization with mean load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NLP for variance optimization
% --->
tic;
ineq_A = dll';
ineq_b = param.vol;
vec_lb = sparse(nm,1);
vec_ub = Inf * ones(nm,1);
x_k = nom_x + (5* ones(nm,1));

options = optimoptions('fmincon','Display','final-detailed',...
    'MaxFunctionEvaluations', 10^10, 'MaxIterations', 10^5,...
    'OptimalityTolerance', 10^(-6), 'Algorithm','interior-point');

[x_k, f_val, exitflag, output]...
    = fmincon(@obj_variance, x_k, ineq_A, ineq_b, [], [], vec_lb, vec_ub, [], options);

var_x = x_k;

elapsed_time = toc;
% <---
% NLP for variance optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing
% --->
fprintf(' ============================================= \n');
fprintf(' ==  Variance optimization >>> \n');
fprintf('     #members = %g \n', nm);
fprintf('     #DOF     = %g \n', nd);
fprintf('     volume = %3.6e [cm^3] <= %3.6e [mm^3] \n',...
    (dll' * var_x * 100), (param.vol * (10^5)) );
opt_matK = sqrtK * diag(var_x) * sqrtK';
[opt_w_worst,~] = worst_analysis(opt_matK, mat_aMM, til_f);
max_pi = til_f' * pinv(opt_matK) * til_f;
fprintf('     compliance for mean force = %3.4f [J] \n',...
    (max_pi * 10) );
fprintf('     worst-case compliance = %3.3f [J] \n',...
    (opt_w_worst * 10) );
fprintf('     discrepancy of compliance = %3.3f [J] \n',...
    (opt_w_worst - max_pi) * 10 );
fprintf('     variance of compliance = %3.3f [J^2] \n',...
    obj_variance(var_x) * 10^2 );
fprintf('     alpha = %3.3e \n',...
    param.alpha );
fprintf('     elapsed time = %4.1f s \n', elapsed_time);
fprintf(' ============================================= \n');
% <---
% Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pareto solutions
% --->
fprintf(' ============================================= \n');
fprintf(' ==  Pareto optimal solutions >>> \n');
his_w_mean  = zeros(1,param.num_Pareto);
his_w_worst = zeros(1,param.num_Pareto);
his_w_discr = zeros(1,param.num_Pareto);
his_w_appVar = zeros(1,param.num_Pareto);
his_w_quadVar = zeros(1,param.num_Pareto);
his_x = cell(1,param.num_Pareto);

load('nominal_opt.mat');
x_k = nominal_opt;

delta_pi = (max_pi - min_pi) / (param.num_Pareto-1);
for iP = 1:param.num_Pareto
    bar_pi = min_pi + ((iP-1) * delta_pi);
    
    if iP > 1
        x_k = x_k + (1 * ones(nm,1));
        
        [x_k, f_val, exitflag, output]...
            = fmincon(@obj_variance, x_k, ineq_A, ineq_b, [], [],...
            vec_lb, vec_ub, @cstr_compliance, options);
    end
    
    x_k = max(0, x_k);
    Par_x = x_k;

    Par_matK = sqrtK * diag(Par_x) * sqrtK';
    [Par_w_worst,~] = worst_analysis(Par_matK, mat_aMM, til_f);
    Par_u = pinv(Par_matK) * til_f;
    Par_pi = til_f' * Par_u;
    appVar = 4 * Par_u' * matSigma * Par_u;
    quadVar = appVar + (2 * trace( (matSigma * pinv(Par_matK))^2 ) );
    fprintf('     compliance for mean force = %3.4f [J] <= %3.4f \n',...
        (Par_pi * 10), bar_pi);
    fprintf('     worst-case compliance = %3.3f [J] \n',...
        (Par_w_worst * 10));
    fprintf('     discrepancy of compliance = %3.3f [J] \n',...
        (Par_w_worst - Par_pi) * 10 );
    fprintf('     variance of compliance = %3.3f [J^2] \n',...
        obj_variance(Par_x) * 10^2 );
    %
    his_w_mean(iP)  = Par_pi * 10;
    his_w_worst(iP) = Par_w_worst * 10;
    his_w_discr(iP) = (Par_w_worst - Par_pi) * 10;
    his_w_appVar(iP) = appVar * 10^2;
    his_w_quadVar(iP) = quadVar * 10^2;
    his_x{iP} = Par_x;
    fprintf(' ---------- \n');
end
fprintf(' ============================================= \n');
figure;
plot(his_w_mean, his_w_worst, 'bo-', 'LineWidth',1.5);
hold on;
xlim([72, 74]);
ylim([88.2, 89.2]);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('Nominal-case compliance (J)', 'Interpreter', 'latex');
ylabel('Worst-case compliance (J)', 'Interpreter', 'latex');
if Flag.save == 1
    saveas(gcf, strcat('varia_Pareto_1_alpha',num2str(param.alpha)), 'epsc');
end



figure;
plot(his_w_mean, his_w_discr, 'rx-', 'LineWidth',1.5);
hold on;
xlim([72, 74]);
ylim([15.4, 16.6]);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('Nominal-case compliance (J)', 'Interpreter', 'latex');
ylabel('Discrepancy of compliance (J)', 'Interpreter', 'latex');
if Flag.save == 1
    saveas(gcf, strcat('varia_Pareto_2_alpha',num2str(param.alpha)), 'epsc');
end


figure;
plot(his_w_appVar, his_w_discr, ':', 'Color',[0.7 0.4 0.3], 'LineWidth',2);
hold on;
plot(his_w_quadVar, his_w_discr, 'g-', 'LineWidth',1.5);
% xlim([23, 25.5]);
ylim([15.4, 16.6]);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('Variance of compliance (J$^2$)', 'Interpreter', 'latex');
ylabel('Discrepancy of compliance (J)', 'Interpreter', 'latex');
if Flag.save == 1
    saveas(gcf, strcat('varia_Pareto_3_alpha',num2str(param.alpha)), 'epsc');
end

figure;
plot(his_w_mean, his_w_quadVar, 'gs-', 'LineWidth',1.5);
hold on;
xlim([72, 74]);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('Nominal-case compliance (J)', 'Interpreter', 'latex');
ylabel('Variance of compliance (J$^2$)', 'Interpreter', 'latex');
if Flag.save == 1
    saveas(gcf, strcat('varia_Pareto_4_alpha',num2str(param.alpha)), 'epsc');
end


[~] = draw_cs_specified_width(coord_x, irr, his_x{end});
if Flag.save == 1
    saveas(gcf, strcat(filename.std,'_varia_alpha',num2str(param.alpha)), 'epsc');
end
% <---
% Pareto solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

