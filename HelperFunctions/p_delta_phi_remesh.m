% Extract phase difference and amplitude difference for a given point
function [d_phi,d_peak] = p_delta_phi_remesh(point)
  
    old_mesh = length(point.mesh);
    new_mesh = old_mesh*10; % refine mesh by a factor of 10 - probably not needed
  
    point = p_remesh(point,3,new_mesh);
    t = point.mesh;
    X = point.profile;
    vec = {X(1,:),X(3,:)};
    peak = zeros(1,2);
    amp = zeros(1,2);

    x_inds = [1,3];
  
    for i = 1:2
      x = vec{i};
      f_p = @(j) fourier_deriv(t,x,j);
      [vals,loc] = findpeaks(x);
      if length(vals) > 1
        [~,id] = max(vals);
        loc = loc(id);
      elseif isempty(vals)
        [~,loc] = max(x);
      end
      tp = t(loc);
      try
        peak(i) = fzero(f_p,round(tp,2));
        amp(i) = interp1(t,X(x_inds(i),:), peak(i));
      catch
        peak(i) = 0 ;
        amp(i) = interp1(t,X(x_inds(i),:), peak(i));
      end
    end
  
    d_peak = abs(amp(2) - amp(1))/min(amp);
  
    diff = peak(2)-peak(1);
    diff = mod(diff,1);
    if diff > 0.5
      diff = 1-diff;
    end
    d_phi = diff;

end