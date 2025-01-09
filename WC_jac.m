function J = WC_jac(y, p)
  % Jacobian of Wilson-Cowan model
  J = zeros(2,2);
  E = y(1);
  I = y(2);

  f = @(u) 1.0/(1.0+exp(-u));
  fp = @(u) f(u).*(1.0-f(u));

  J(1,1) = -1.0 + p.a*fp(p.p_e+p.a*E-p.b*I);
  J(1,2) = -p.b*fp(p.p_e+p.a*E-p.b*I);
  J(2,1) = -p.c*fp(p.p_i+p.c*E-p.d*I);
  J(2,2) = -1.0 - p.d*fp(p.p_i+p.c*E-p.d*I);
end