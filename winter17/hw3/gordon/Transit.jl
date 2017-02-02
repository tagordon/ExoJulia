module Transit
export transit

function transit(Rp, Rstar, R)
  trans_curve = zeros(length(R))
  for (i, r) in enumerate(R)
    if r > (Rstar + Rp)
      Reff = 0
    elseif r < (Rstar - Rp)
      Reff = Rp
    else
      x = ((r^2) - (Rp^2 - Rstar^2))/(2*r)
      cos_thetap = (r - x)/Rp
      cos_thetastar = x/Rstar
      thetap = acos(cos_thetap)
      thetastar = acos(cos_thetastar)
      Ap = (Rp^2)*thetap - ((Rp^2)/2)*sin(2*thetap)
      Astar = (Rstar^2)*thetastar - ((Rstar^2)/2)*sin(2*thetastar)
      A = Ap + Astar
      Reff = sqrt(A/(pi))
    end
    trans_curve[i] = 1-(Reff/Rstar)^2
  end
  return trans_curve
end
end
