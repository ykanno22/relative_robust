function [comp_cstr, ceq] = cstr_compliance(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2019 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global sqrtK til_f bar_pi

matK = sqrtK * diag(x) * sqrtK';
matK = matK + (10^(-12)*speye(size(matK,1)));

vec_u = matK \ til_f;

comp_cstr = (til_f' * vec_u) - bar_pi;

ceq = [];