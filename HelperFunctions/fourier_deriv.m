% Test Fourier script

function grad = fourier_deriv(t,X,point)
  
  t_discrete = linspace(t(1), t(end), 1024)';
  
  X = interp1(t,X,t_discrete);
  signal = X;
  
  nFCs = 10;
  n = 0:nFCs-1;
  Tp = t_discrete(end)-t_discrete(1);
  
  phases = 2*pi*t_discrete*n/Tp;
  S = sin(phases);
  C = cos(phases);
  
  SIGNAL = repmat(signal, [1,nFCs]);
  
  INTEGRANDsin = SIGNAL.*S;
  INTEGRANDcos = SIGNAL.*C;
  
  FCsin = 2/(1023)*trapz(INTEGRANDsin, 1);
  FCcos = 2/(1023)*trapz(INTEGRANDcos, 1);
  
  % Project back into real space
  
  t_alt = point;
  
  phases_alt = 2*pi*t_alt*n/Tp;
  S_alt = cos(phases_alt);
  C_alt = -sin(phases_alt);
  
  FCsin_deriv = 2*pi*n/Tp.*FCsin;
  FCcos_deriv = 2*pi*n/Tp.*FCcos;
  SFCs_rep = repmat(FCsin_deriv, [length(point),1]);
  CFCs_rep = repmat(FCcos_deriv, [length(point),1]);
  
  grad = sum(SFCs_rep.*S_alt, 2) + sum(CFCs_rep.*C_alt,2);

end
