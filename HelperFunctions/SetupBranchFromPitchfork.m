% Setup branch continuation from pitchfork

function [nfuncs,nbranch,suc] = SetupBranchFromPitchfork(sol_type, funcs, branch, ind_pt, ind_cont)

  switch sol_type
  
    case 'in-phase'

      % Shift for additional constraint
      shift = [0,1];

      % Integral boundaries for additional constraint
      condprojint = linspace(0,1,6)'*[1,1];
  
    case 'anti-phase'

      % Shift for additional constraint
      shift = [1,2];

      % Integral boundaries for additional constraint
      condprojint = linspace(0,0.5,6)'*[1,1];
  
  end
  
  % System size
  xdim = 4; % In our case, this is 2x2 = 4

  % Create branch
  Rsym = zeros(4);
  Rsym(1,3) = 1;
  Rsym(2,4) = 1;
  Rsym(3,1) = 1;
  Rsym(4,2) = 1;

  sbxsym = @(p) dde_psol_lincond(p, xdim, 'x', ...
    'trafo', Rsym, 'shift', shift, 'condprojint', condprojint);
  poev1args = {'usercond', {sbxsym}};

  % Make sure arguments are only passsed on to find pitchfork and not to
  % the actual solutions themselves
  nspoev1args = addprefix('SetupPOEV1', poev1args);

    % end of "this"
  [nfuncs,nbranch,suc] = SetupPsol(funcs, branch, ...
      ind_pt, 'branch_off', 'POEV1',...
      'print_residual_info', 1, 'outputfuncs', true, nspoev1args{:}, ...
      'contpar',ind_cont, ...
      'continuation.plot_measure', ...
          {@(p) p.parameter(ind_cont), @(p) p_delta_phi_remesh(p)});
  
  % Follow branch
  nbranch.parameter.min_bound = [ind_cont, 0.0];
  nbranch.parameter.max_bound = [ind_cont, pi];

end