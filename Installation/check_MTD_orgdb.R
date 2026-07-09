#!/usr/bin/env Rscript

required_packages <- c(
    "AnnotationDbi",
    "AnnotationForge",
    "biomaRt",
    "GenomeInfoDb",
    "GO.db",
    "DBI",
    "RSQLite"
)

package_status <- vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
)

cat("============================================================\n")
cat("MTD_orgdb environment validation\n")
cat("============================================================\n")

cat("R version:\n")
cat("  ", R.version.string, "\n", sep = "")

cat("R home:\n")
cat("  ", R.home(), "\n", sep = "")

cat("R library paths:\n")
for (library_path in .libPaths()) {
    cat("  ", library_path, "\n", sep = "")
}

cat("\nPackage status:\n")

for (package_name in names(package_status)) {
    installed <- package_status[[package_name]]

    version <- if (installed) {
        as.character(utils::packageVersion(package_name))
    } else {
        "not installed"
    }

    result <- if (installed) "PASS" else "FAIL"

    cat(
        sprintf(
            "  %-20s %-4s  %s\n",
            package_name,
            result,
            version
        )
    )
}

if (!all(package_status)) {
    missing_packages <- names(package_status)[!package_status]

    stop(
        "MTD_orgdb is missing required packages: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

cat("\n[OK] MTD_orgdb is ready for OrgDb construction.\n")
