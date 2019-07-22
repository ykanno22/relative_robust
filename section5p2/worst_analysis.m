function [post_w, post_y] = worst_analysis(matK, matbUSU, til_f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2019 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cvx_begin sdp quiet
    cvx_solver sdpt3
    cvx_precision best
    variable post_y(1,1)
    variable post_w(1,1)
    minimize( post_w )
    subject to
    [post_w, 1, til_f';
        1, post_y, sparse(1,size(matK,1));
        til_f, sparse(size(matK,1),1), matK - (post_y * matbUSU)] >= 0;
cvx_end
