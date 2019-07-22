%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2019 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear;
close all;
%
Flag.save = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
param.nominal_load = 10.0;
param.alpha = 1.0;
param.vol = 100; % *(10^5) [mm^3]
param.num_Pareto = 30;
%
eeee = 2.0;
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data of the structure
% --->
nx = 3;
ny = 3;

filename.std = strcat('x',num2str(nx),'_y',num2str(ny));

[dll,matH,coord_x,ir,irr,ird] = member(nx, ny, Inf);

nk = size(coord_x,1); num.node   = nk;
nd = size(matH,1);    num.degree = nd;
nm = size(matH,2);    num.member = nm;

[~] = draw_cs_ini(coord_x, irr, ones(nm,1), nx, ny);
% if Flag.save == 1
%     saveas(gcf, strcat(filename.std,'_initial'), 'epsc');
% end
% <---
% Data of the structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vector
% --->
Idx.deg_of_load = 2 * nx;

til_f = sparse(zeros(nd,1));
til_f(Idx.deg_of_load) = -param.nominal_load;
%
matE = [1, -1; -1, 0.3];
matE = matE / max(abs(eig(matE)));
matM = zeros(nd,2);
matM((Idx.deg_of_load-1):(Idx.deg_of_load),:) = matE;
matM = sparse(matM);
matSigma = ((param.alpha/3)^2) * matM * (matM');
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

[~] = draw_cs_specified_width(coord_x, irr, nom_x);
if Flag.save == 1
    saveas(gcf, strcat(filename.std,'_nominal_mean'), 'epsc');
end
fprintf(' ============================================= \n');
% <---
% Nominal optimization with mean load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CCP for discrepancy optimization
% --->
x_k = ones(nm,1);
diff_x = Inf;
Iter = 0;

matK_k = sqrtK * diag(x_k) * sqrtK';
u_k = matK_k \ til_f;
grad_pi_k = -(sqrtK' * u_k) .* (sqrtK' * u_k);

tic;
while (diff_x > 10^(-3)) && (Iter < 5000)
    cvx_begin sdp quiet
        cvx_solver sedumi
        cvx_precision best
        variable rel_x(nm,1)
        variable rel_y(1,1)
        variable rel_z(1,1)
        minimize( rel_z - (grad_pi_k' * rel_x) +  10^(-4)*norm(rel_x - x_k) )
        subject to
            [rel_z, 1, til_f';
                1, rel_y, zeros(1,nd);
                til_f, zeros(nd,1),...
                (sqrtK * diag(rel_x) * sqrtK') - (rel_y * mat_aMM)] >= 0;
            dll' * rel_x <= param.vol;
            rel_x >= 0;
    cvx_end
    
    diff_x = norm(x_k - rel_x);
    x_k = rel_x;
    Iter = Iter + 1;
    
    matK_k = sqrtK * diag(x_k) * sqrtK';
    u_k = matK_k \ til_f;
    grad_pi_k = -(sqrtK' * u_k) .* (sqrtK' * u_k);
    [rel_w_worst,~] = worst_analysis(matK_k, mat_aMM, til_f);
    fprintf('      Iter %g: obj=%8.4f  nom.comp.=%8.4f change=%7.3f\n',...
        Iter, rel_w_worst - (til_f'*u_k), til_f'*u_k, diff_x);
end
elapsed_time = toc;
% <---
% CCP for discrepancy optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing
% --->
fprintf(' ============================================= \n');
fprintf(' ==  Discrepancy optimization >>> \n');
fprintf('     #members = %g \n', nm);
fprintf('     #DOF     = %g \n', nd);
fprintf('     volume = %3.6e [cm^3] <= %3.6e [mm^3] \n',...
    (dll' * rel_x * 100), (param.vol * (10^5)) );
rel_matK = sqrtK * diag(rel_x) * sqrtK';
[rel_w_worst,~] = worst_analysis(rel_matK, mat_aMM, til_f);
max_pi = til_f' * (rel_matK \ til_f);
fprintf('     compliance for mean force = %3.4f [J] \n',...
    (max_pi * 10) );
fprintf('     worst-case compliance = %3.3f [J] ; obj. = %3.3f \n',...
    (rel_w_worst * 10), rel_z );
fprintf('     discrepancy of compliance = %3.3f [J] \n',...
    (rel_w_worst - max_pi) * 10 );
fprintf('     alpha = %3.3e \n',...
    param.alpha );
fprintf('     elapsed time = %4.1f s \n', elapsed_time);
fprintf(' ============================================= \n');

[~] = draw_cs_specified_width(coord_x, irr, rel_x);
if Flag.save == 1
    saveas(gcf, strcat(filename.std,'_discr_alpha',num2str(param.alpha)), 'epsc');
end
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
% 
x_k = nom_x;
%
delta_pi = (max_pi - min_pi) / (param.num_Pareto-1);
for iP = 1:param.num_Pareto
    bar_pi = min_pi + ((iP-1) * delta_pi);

    diff_x = Inf;
    x_k = x_k + (5* ones(nm,1));
    Iter = 0;
    
    matK_k = sqrtK * diag(x_k) * sqrtK';
    u_k = matK_k \ til_f;
    grad_pi_k = -(sqrtK' * u_k) .* (sqrtK' * u_k);

    while (diff_x > 10^(-3)) && (Iter < 5000)
        cvx_begin sdp quiet
            cvx_solver sedumi
            cvx_precision best
            variable Par_x(nm,1)
            variable Par_y(1,1)
            variable Par_z(1,1)
            minimize( Par_z - (grad_pi_k' * Par_x) +  10^(-4)*norm(Par_x - x_k) )
            subject to
                [Par_z, 1, til_f';
                    1, Par_y, zeros(1,nd);
                    til_f, zeros(nd,1),...
                    (sqrtK * diag(Par_x) * sqrtK') - (Par_y * mat_aMM)] >= 0;
                [bar_pi, til_f';
                    til_f, (sqrtK * diag(Par_x) * sqrtK')] >= 0;
                dll' * Par_x <= param.vol;
                Par_x >= 0;
        cvx_end
    
        diff_x = norm(x_k - Par_x);
        x_k = Par_x;
        Iter = Iter + 1;
        
        matK_k = sqrtK * diag(x_k) * sqrtK';
        u_k = matK_k \ til_f;
        grad_pi_k = -(sqrtK' * u_k) .* (sqrtK' * u_k);
        [Par_w_worst,~] = worst_analysis(matK_k, mat_aMM, til_f);
        fprintf('      Iter %g: obj=%8.4f  nom.comp.=%8.4f<=%8.4f change=%7.3f\n',...
            Iter, Par_w_worst - (til_f'*u_k), til_f'*u_k, bar_pi, diff_x);
    end

    Par_matK = sqrtK * diag(Par_x) * sqrtK';
    [Par_w_worst,~] = worst_analysis(Par_matK, mat_aMM, til_f);
    Par_u = Par_matK \ til_f;
    Par_pi = til_f' * Par_u;
    appVar = 4 * Par_u' * matSigma * Par_u;
    quadVar = appVar + (2 * trace( (matSigma * pinv(Par_matK))^2 ));
    fprintf('     compliance for mean force = %3.4f [J] <= %3.4f \n',...
        (Par_pi * 10), bar_pi);
    fprintf('     worst-case compliance = %3.3f [J] ; obj. = %3.3f \n',...
        (Par_w_worst * 10), Par_z );
    fprintf('     discrepancy of compliance = %3.3f [J] \n',...
        (Par_w_worst - Par_pi) * 10 );
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
plot(his_w_mean, his_w_worst, 'b-', 'LineWidth',1.5);
hold on;
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
% xlim([360, 400]);
% ylim([510, 535]);
ylim([389, 390.2]);
xlabel('Nominal-case compliance (J)', 'Interpreter', 'latex');
ylabel('Worst-case compliance (J)', 'Interpreter', 'latex');
if Flag.save == 1
    saveas(gcf, strcat(filename.std,'_discr_Pareto_1_alpha',num2str(param.alpha)), 'epsc');
end

figure;
plot(his_w_mean, his_w_discr, 'r-', 'LineWidth',1.5);
hold on;
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('Nominal-case compliance (J)', 'Interpreter', 'latex');
ylabel('Discrepancy of compliance (J)', 'Interpreter', 'latex');
if Flag.save == 1
    saveas(gcf, strcat(filename.std,'_discr_Pareto_2_alpha',num2str(param.alpha)), 'epsc');
end



figure;
plot(his_w_appVar, his_w_discr, ':', 'Color',[0.7 0.4 0.3], 'LineWidth',2);
hold on;
plot(his_w_quadVar, his_w_discr, 'g-', 'LineWidth',1.5);
xlim([45, 85]);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlabel('Variance of compliance (J$^2$)', 'Interpreter', 'latex');
ylabel('Discrepancy of compliance (J)', 'Interpreter', 'latex');
if Flag.save == 1
    saveas(gcf, strcat(filename.std,'_discr_Pareto_3_alpha',num2str(param.alpha)), 'epsc');
end



[~] = draw_cs_specified_width(coord_x, irr, his_x{1});
% <---
% Pareto solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

