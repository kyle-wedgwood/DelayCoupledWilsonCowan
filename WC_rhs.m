function dydt = WC_rhs(~, y, vp, fp)

  % Unpack variables
  E = y(1);
  I = y(2);

  % RHS for Wilson-Cowan model
  f = @(u) 1.0./(1.0+exp(-u));

  dydt = [ -E + f(vp(1)+fp.a*E-fp.b*I);
           -I + f(vp(2)+fp.c*E-fp.d*I)]; 

end