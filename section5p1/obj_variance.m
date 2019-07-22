function obj_val = obj_variance(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2019 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global sqrtK til_f matSigma

matK = sqrtK * diag(x) * sqrtK';
matK = matK + (10^(-12)*speye(size(matK,1)));

vec_u = matK \ til_f;
appVar = 4 * vec_u' * matSigma * vec_u;
obj_val = appVar + (2 * trace( (matSigma * pinv(matK))^2 ) );

