# Calculate inhibition
# Matrix approach

inhib <- function(im0, 
                  im1, 
                  iconc, 
                  combine = prod) {

  # Extract groups and inhibitors for dimensions but mainly for debugging
  inhibs <- colnames(im0)
  grps <- rownames(im0)

  # Convert concentrations to matrix
  cmat <- matrix(rep(iconc, length(grps)), 
                 nrow = length(grps), 
                 byrow = TRUE,
                 dimnames = list(grps, 
                                 inhibs))
  # Slope matrix
  smat <- (1 - 0) / (im1 - im0)

  # Reduction matrix
  rmat <- 1 - (cmat - im0) * smat
  rmat[rmat > 1] <- 1
  rmat[rmat < 0] <- 0

  # Get product by microbial group
  red <- apply(rmat, 1, combine)

}
