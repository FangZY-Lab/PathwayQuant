# Declare non-standard-evaluation symbols used by limma::makeContrasts and
# base::subset so that R CMD check does not flag them as undefined globals.
utils::globalVariables(c("H", "L", "logFC"))
