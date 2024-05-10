# ABM
R model for anaerobic microbial degradation (biodegradation) of organic matter with multiple microbial groups

# Description
The ABM model predicts conversion of animal manure or other high-moisture organic wastes to methane and carbon dioxide under anaerobic conditions.
The name comes from **a**naerobic **b**iodegradation **m**odel. 
With multiple microbial groups and group-specific parameters describing kinetics and yield, the model can predict realistic short- and long-term responses to temperature change.

# Maintainers
Frederik Dalby (<https://au.dk/fd@bce.au.dk>) and Sasha D. Hafner (<https://au.dk/sasha.hafner@bce.au.dk>).

# Installation and getting started
This R package can be installed using `devtools::install_github("AU-BCE-EE/ABM", build_vignettes = TRUE)`.
Load the package with `library(ABM)`, and open up the vignette with `vignette("ABM_start")`.

If you have trouble installing ABM or devtools, make sure that you are working with the latest version of R. If you have updated R but still have trouble, please contact us or submit a new issue, explaining the problem.

# Problems?
Report problems or ask questions using the Issues page.

# Research papers
See these papers for more information on the model, as well as examples of applications.

Dalby, F.R., Hansen, M.J., Guldberg, L.B., Hafner, S.D., Feilberg, A., 2023. Simple Management Changes Drastically Reduce Pig House Methane Emission in Combined Experimental and Modeling Study. Environ. Sci. Technol. <https://doi.org/10.1021/acs.est.2c08891> ([download here](https://drive.google.com/file/d/1FsFzUrbM2O3GlTm3TYHjrlTfr_mFfLRC/view?usp=share_link))

Dalby, F.R., Hansen, M.J., Guldberg, L.B., Hafner, S.D., Feilberg, A., 2023. Simple management changes drastically reduce pig house methane emission in combined experimental and modeling study. Environ. Sci. Technol. <https://doi.org/10.1021/acs.est.2c08891> ([download here](https://drive.google.com/file/d/1zTHmyg3EHHqtb8RMoA6VVf_DnNxs-y-k/view?usp=drive_link)).

