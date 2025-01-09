function dudt = WC_rhs_coupled(xx, pars, fp)
  % Adding coupling to model and adjusting for use with dde-biftool

  % Extract parameters
  gain = pars(1);

  % Extract variables at current time
  EI_1 = xx([1,2],1);
  EI_2 = xx([3,4],1);

  % Extract lagged variables (only voltages)
  Elag_1 = xx(1,2);
  Elag_2 = xx(3,2);

  % Evaluate coupling term
  coupling = gain*[Elag_2-EI_1(1); 0.0; Elag_1-EI_2(1); 0.0];

  % Evaluate uncoupled model
  vp = [fp.p_e; fp.p_i];
  dudt = [ WC_rhs(0, EI_1, vp, fp);
           WC_rhs(0, EI_2, vp, fp)];

  % Add in coupling
  dudt = dudt + coupling;

end