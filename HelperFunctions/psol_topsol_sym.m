function [funcs,branch,suc] = psol_topsol_sym(funcs,p0,orbit,Tp,col_degree,nr_int,ind_cont,step,sym_cond,method)
%% create starting guess for periodic solution derived from periodic solution with additional constraints
% function [psol_point,stpcond]=psol_topsol(funcs,point,ampl,col_degree,nr_int)
% INPUT:
%   funcs problem functions
%   p0 initial parameter vector
%   orbit interpolant for periodic orbit to start frm
%   Tp period of orbit
%	col_degree piecewise polynomial degree for periodic solution
%	nr_int number of intervals for mesh
%   ind_cont paramter to continue in

% OUTPUT:
%	psol_point starting guess for periodic solution derived from point
%	stpcond steplength condition for use in correction of guess
% 
%
    % if no extra condition, do not set flag
    branch = df_brnch(funcs, ind_cont, 'psol');

    mesh = linspace(0, 1, nr_int*col_degree+1);
    profile = orbit(mesh*Tp)';

    point.kind = 'psol';
    point.parameter = p0;
    point.mesh = mesh;
    point.degree = col_degree;
    point.profile = profile;
    point.period = Tp;
    
    if exist('method', 'var')
      branch.method = method;
    end

    % Apply extra condition here
    funcs = dde_funcs_add_cond(funcs, sym_cond);
    branch.method.point.extracolumns = 'auto';

    % First point
    [point,suc] = p_correc(funcs, point, [], [], branch.method.point, 0);
    point.stability = p_stabil(funcs, point, branch.method.stability);

    branch = rmfield(branch, 'point');
    branch.point(1) = point;

    % Second point
    point2 = point;
    %point2.parameter = p0 + [0.1*(ind_cont==indtau), 0.01*(ind_cont==indeps), 0.01*(ind_cont==3)]; %adjust for general params
    % One percent perturbation in direction of cont par
    pd = zeros(size(p0));

    if ~exist('step', 'var')
      step = 0.01;
    end
    if p0(ind_cont) == 0
        pd(ind_cont) = step;
    else
        pd(ind_cont) = step*abs(p0(ind_cont));
    end
    point2.parameter = p0 + pd;

%     [point2,suc] = p_correc(funcs, point2, ind_cont, [], branch.method.point, 0, point);
    [point2,suc] = p_correc(funcs, point2, [], [], branch.method.point, 0, point);

    % compute stability at point
    point2.stability = p_stabil(funcs, point2, branch.method.stability);
    branch.point(2) = point2;

end