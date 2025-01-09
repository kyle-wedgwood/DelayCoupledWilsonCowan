function [funcs,branch,suc] = SetupBranchSym(sol_type, funcs, pars, ind_cont, mesh, profile, Tp)

  switch sol_type
  
    case 'in-phase'
  
      % Copy orbit - start from in-phase solution
      coupled_orbit = @(s) [interp1(mesh, profile, s), ...
        interp1(mesh, profile, s)];
  
      % Discretisation options - in-phase solution
      col_deg  = 3;
      nr_int   = 30;
  
      % Step size for initial plotting - in-phase solution
      step = 0.001;

      % Shift for additional constraint
      shift = [0,1];

      % Integral boundaries for additional constraint
      condprojint = linspace(0,1,6)'*[1,1];
  
    case 'anti-phase'
      % Copy orbit - start from antiphase solution
      coupled_orbit = @(s) [interp1(mesh, profile, s), ...
        interp1(mesh, profile, s-0.5*Tp+(s<0.5*Tp)*Tp)];
  
      % Discretisation options - antiphase solution
      col_deg  = 4;
      nr_int   = 30;
  
      % Step size for initial plotting - antiphase solution
      step = 0.001;

      % Shift for additional constraint
      shift = [1,2];

      % Integral boundaries for additional constraint
      condprojint = linspace(0,0.5,6)'*[1,1];

    case 'out-of-phase'
      % Copy orbit - start from antiphase solution
      coupled_orbit = @(s) [interp1(mesh, profile, s), ...
        interp1(mesh, profile, s-0.1*Tp+(s<0.1*Tp)*Tp)];
  
      % Discretisation options - antiphase solution
      col_deg  = 4;
      nr_int   = 30;
  
      % Step size for initial plotting - antiphase solution
      step = 0.001;
  
  end
  
  % Create branch
  Rsym = zeros(4);
  Rsym(1,3) = 1;
  Rsym(2,4) = 1;
  Rsym(3,1) = 1;
  Rsym(4,2) = 1;
 
  xdim = 4; % size of system. In my case this is 2x2 = 4
  psolsym = @(p) dde_psol_lincond(p, xdim, 'profile', 'trafo', Rsym, ...
      'shift', shift,...
      'condprojint', condprojint);
  [funcs,branch,suc] = psol_topsol_sym(funcs, pars, coupled_orbit, ...
    Tp, col_deg, nr_int, ind_cont, step, psolsym);

  % Plotting options
  branch.method.continuation.plot_measure = ...
    {@(p) p.parameter(ind_cont), @(p) p_delta_phi_remesh(p)};
  
  % Follow branch
  branch.parameter.min_bound = [ind_cont, 0.0];
  branch.parameter.max_bound = [ind_cont, pi];

end