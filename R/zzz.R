.onLoad <- function(libname, pkgname) {
  if (exists(".on_load_options", envir = asNamespace(pkgname))) {
    get(".on_load_options", envir = asNamespace(pkgname))(libname, pkgname)
  }

  methods <- list(
    print = print.parametric_hrf_fit,
    coef = coef.parametric_hrf_fit,
    summary = summary.parametric_hrf_fit,
    fitted = fitted.parametric_hrf_fit,
    residuals = residuals.parametric_hrf_fit
  )
  for (generic in names(methods)) {
    utils::registerS3method(generic, "parametric_hrf_fit", methods[[generic]],
                            envir = parent.env(environment()))
  }
  utils::registerS3method("print", "summary_parametric_hrf_fit",
                          print.summary_parametric_hrf_fit,
                          envir = parent.env(environment()))
}
