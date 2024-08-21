\name{parameters}

\alias{parameters}
\alias{ABM_parameters}
\alias{arrh_pars2.0}
\alias{chem_pars2.0}
\alias{grp_pars2.0}
\alias{man_pars2.0}
\alias{mic_pars2.0}
\alias{wthr_pars2.0}
\alias{arrh_pars1.0}
\alias{chem_pars1.0}
\alias{grp_pars1.0}
\alias{man_pars1.0}

\docType{data}

\title{
  Default parameter sets for ABM
}

\description{
  Default parameter sets for the \code{\link{abm}} function.
}

\usage{

  data(grp_pars2.0)
  data(grp_pars1.0)
  data(man_pars2.0)
  data(man_pars1.0)
  data(mic_pars2.0)
  data(wthr_pars2.0)
  data(chem_pars2.0)
  data(arrh_pars2.0)
  
}

\format{

  These are default parameter sets.
  Sets with \code{2.0} are the latest default sets.
  The \code{1.0} versions are from Dalby et al. (2023b) (see Source section). 
  The \code{arrh_pars2.0} set contains hydrolysis rates of CP, RFd and starch and sugars, and those constants come from Dalby et al. (2023b).
  The values in \code{man_pars2.0} are also derived from the measurements Dalby et al. (2023b) but are not presented in that paper. 
  The \code{grp_pars2.0} set contains methanogen substrate uptake rates (qhat_opt), 
  which come from fitting to multiple full scale datasets from pig farms. This data is not yet published.
  The \code{grp_pars2.0} also contains inhibition kinetics, which is described Astals et al. (2018) for NH3 and NH4. Upper inhibition threshold for NH4 (ki_NH4_max) has been changed compared to Astals et al. (2018) by fitting to multiple full scale datasets from pig farms. H2S inhibition constants are derived from data in O'Flaherty et al. (1998). 
  
  
}

\source{
Dalby, F.R., Hafner, S.D., Petersen, S.O., Vanderzaag, A., Habtewold, J., Dunfield, K., Chantigny, M.H., Sommer, S.G., 2021. A mechanistic model of methane emission from animal slurry with a focus on microbial groups. PLOS ONE 16, e0252881. \doi{10.1371/journal.pone.0252881}

Dalby, F.R., Ambrose, H.W.,  Poulsen, J.S., Nielsen, J.L., Adamsen, A.P.S., 2023a. Pig slurry organic matter transformation and methanogenesis at ambient storage temperatures. JEQ. \doi{10.1002/jeq2.20512}

Dalby, F.R., Hansen, M.J., Guldberg, L.B., Hafner, S.D., Feilberg, A., 2023b. Simple Management Changes Drastically Reduce Pig House Methane Emission in Combined Experimental and Modeling Study. Environ. Sci. Technol. \doi{10.1021/acs.est.2c08891}

Astals, S., Peces, M., Batstone, D. J., Jensen, P. D., Tait, S. 2018. Characterising and modelling free ammonia and ammonium inhibition in anaerobic systems. Water Research. \doi{org/10.1016/j.watres.2018.06.021}

O'Flaherty, V., Mahony T., O'Kennedy, R., Colleran, E. 1998. Effect of pH on growth kinetics and sulphide toxicity thresholds of a range of methanogenic syntrophic and sulphate-reducing bacteria. Vol. 33, 5. p555-569. 

}

\examples{
  out1 <- abm(man_pars = man_pars2.0, grp_pars = grp_pars1.0)
}

\keyword{datasets}
