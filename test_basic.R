library(devtools)
load_all()
data_obj <- matrix(rnorm(20), nrow = 10, ncol = 2)
event_obj <- matrix(rbinom(10, 1, 0.2), ncol = 1)
res <- estimate_parametric_hrf(data_obj, event_obj, global_refinement = FALSE, verbose = FALSE)
cat("SUCCESS: Function worked\!\n")
cat("Result class:", class(res), "\n")
EOF < /dev/null