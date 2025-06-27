# load nothing extra — this only uses base R

# Helpers
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)
lnM_core <- function(Delta, MSW, n0) 0.5 * (log(Delta) - log(n0) - log(MSW))

# Debug Δ‐method: returns raw vs capped variance
lnM_delta_ind_dbg <- function(x1, x2, s1, s2, n1, n2, maxVar = 20) {
  h     <- n1 * n2 / (n1 + n2)
  MSB   <- h * (x1 - x2)^2
  MSW   <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(var_raw = NA_real_, var_capped = NA_real_))
  # compute gradients
  sD2       <- s1^2/n1 + s2^2/n2
  dif       <- x1 - x2
  vB        <- h^2 * (2*sD2^2 + 4*sD2*dif^2)
  vW        <- 2 * MSW^2 / (n1 + n2 - 2)
  g1        <- 0.5 / Delta
  g2        <- -0.5 * MSB / (Delta * MSW)
  var_raw    <- posify(g1^2 * vB + g2^2 * vW)
  var_capped <- pmin(var_raw, maxVar)
  c(var_raw = var_raw, var_capped = var_capped)
}

# A little function to build the “near‐zero‐Δ” means
make_means <- function(n1, n2, eps) {
  MSW  <- 1      # assume s1 = s2 = 1 → pooled MSW = 1
  h    <- n1*n2/(n1+n2)
  # solve h*(x_diff)^2 - MSW = eps  →  x_diff = sqrt((MSW + eps)/h)
  xdiff <- sqrt((MSW + eps)/h)
  # split symmetrically so mean difference = xdiff
  list(x1bar =  xdiff/2,
       x2bar = -xdiff/2,
       s1 = 1, s2 = 1)
}

# Now three examples on (n1, n2) = (3,7) with eps = 1e-4, 1e-6, 1e-8
for (eps in c(1e-4, 1e-6, 1e-8)) {
  mm <- make_means(n1 = 3, n2 = 7, eps = eps)
  res <- lnM_delta_ind_dbg(mm$x1bar, mm$x2bar, mm$s1, mm$s2,
                           n1 = 3, n2 = 7, maxVar = 20)
  cat(sprintf("eps = %-8g → var_raw = % .3e   var_capped = % .3f\n",
              eps, res["var_raw"], res["var_capped"]))
}