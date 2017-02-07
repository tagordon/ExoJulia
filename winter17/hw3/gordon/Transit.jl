module Transit
export transit, depth, ld_transit

# compute the depth of the transit at a given planet location R
function depth(Rp, Rstar, R)
  if Rp > Rstar
    temp = Rp
    Rp = Rstar
    Rstar = temp
  end
  if R > (Rstar + Rp)
    Reff = 0
  elseif R < (Rstar - Rp)
    Reff = Rp
  else
    x = ((R^2) - (Rp^2 - Rstar^2))/(2*R)
    cos_thetap = (R - x)/Rp
    cos_thetastar = x/Rstar
    thetap = acos(cos_thetap)
    thetastar = acos(cos_thetastar)
    Ap = (Rp^2)*thetap - ((Rp^2)/2)*sin(2*thetap)
    Astar = (Rstar^2)*thetastar - ((Rstar^2)/2)*sin(2*thetastar)
    A = Ap + Astar
    Reff = sqrt(A/(pi))
  end
  depth = (Reff/Rstar)^2
  return depth
end

# compute the full transit for the given parameters
# Rp = planet radius
# Rstar = stellar radius
# t0 = tI (first contact)
# dur = transit duration (TIV - TI)
# b = impact parameter
function transit(t, Rp, Rstar, t0, dur, b)
  lc = zeros(length(t))
  if b >= (Rstar + Rp)
    path_length = 0
  else
    path_length = 2*sqrt((Rstar+Rp)^2 - (b)^2)
  end
  for (i,time) in enumerate(t)
    if time > t0
      x = path_length*(time-t0)/dur - path_length/2.
      R = sqrt(x^2 + (b)^2)
      lc[i] = 1 - depth(Rp, Rstar, R)
    else
      lc[i] = 1
    end
  end
  return lc
end

# compute a limb darkened transit
# ld = limb darkening profile (a function specifying intensity vs. mu)
# ld_args = arguments for the limb darkening profile
function ld_transit(t, Rp, t0, dur, b, ld, ld_args)
  # number of layers
  nlayers = 1000

  # compute the radii of each layer
  mu = linspace(0,1,10000)
  Imin = ld(mu[1], ld_args)
  Imax = ld(mu[end], ld_args)
  dI = (Imax - Imin)/nlayers
  Radii = zeros(nlayers)
  # the first layer should have the radius of the star, which is set to 1
  Radii[1] = 1
  i = 1
  # if the profile is flat, all radii are 1
  if Imin == Imax
    Radii = ones(length(Radii))
    dI = 1./nlayers
  else
    # figure out what the radius of each layer should be
    for m in mu
      if ld(m, ld_args) > (Imin + i*dI)
          Radii[i + 1] = 1 - m
          i = i + 1
      end
    end
  end
  # clip off any layers that aren't needed so that we don't have a stack
  # of points in the middle of the star.
  Radii = Radii[1:i]

  # the first layer should actually be a stack of layers deep enough to
  # give an intensity Imin.
  trans = (1 - transit(t, Rp, 1, t0 , dur, b))*Imin
  # compute the transit for each layer and add them together
  for R in Radii
    # an approximation of the velocity of the planet accross the face of the star.
    v = 2*sqrt((1 + b)^2)/dur
    # get the depth from the transit function, then do some stuff
    # to fix it.
    depth = (1 - transit(t, Rp, R, t0 + (1-R)/v, dur - 2*(1-R)/v, b))
    depth = depth*(R^2)*dI
    trans = trans + depth
  end
  return 1 - trans
end


end
