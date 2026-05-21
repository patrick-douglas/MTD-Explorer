#!/usr/bin/env bash
# MTD R package installer - optimized v4
# Target: conda env MTD, R 4.0.3, Bioconductor 3.12
#
# Main idea:
#   - Do NOT mix current CRAN with R 4.0.3.
#   - Use CRAN snapshot 2022-05-15 because tidyverse needs dplyr >= 1.0.9.
#   - Install critical packages in fresh R sessions to avoid namespace cache issues
#     such as "rlang 0.4.11 is already loaded".
#   - Keep UCSC.utils and rgeos optional/outside required validation for this old stack.
#
# Recommended:
#   conda activate MTD
#   CONDA_ENSURE=0 bash update_fix/Install.R.packages.MTD.sh

set -Eeo pipefail

ENV_NAME="${ENV_NAME:-MTD}"
CONDA_ENSURE="${CONDA_ENSURE:-0}"
INSTALL_OPTIONAL="${INSTALL_OPTIONAL:-1}"
REPAIR_BIOC_CORE="${REPAIR_BIOC_CORE:-1}"
REPAIR_TIDYVERSE="${REPAIR_TIDYVERSE:-1}"
CRAN_REPO="${CRAN_REPO:-https://packagemanager.posit.co/cran/2022-05-15}"
CONDA_SH="${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -d "$SCRIPT_DIR/../update_fix/pvr_pkg" ]]; then
  MTD_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
elif [[ -d "$PWD/update_fix/pvr_pkg" ]]; then
  MTD_DIR="$PWD"
elif [[ -n "${dir:-}" && -d "${dir}/update_fix/pvr_pkg" ]]; then
  MTD_DIR="$dir"
else
  MTD_DIR="$HOME/MTD"
fi

PATCH_DIR="${PATCH_DIR:-$MTD_DIR/update_fix/pvr_pkg}"
LOG_DIR="${LOG_DIR:-$MTD_DIR/update_fix/MTD_R_install_logs}"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/Install.R.packages.MTD_optimized_v4.$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG_FILE") 2>&1

log() {
  printf '\n[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
  printf '\n[%s] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2
  exit 1
}

cleanup_locks() {
  if [[ -n "${CONDA_PREFIX:-}" && -d "$CONDA_PREFIX/lib/R/library" ]]; then
    rm -rf "$CONDA_PREFIX/lib/R/library/00LOCK"* || true
  fi
}

run_r() {
  local label="$1"
  shift
  log "$label"
  Rscript --vanilla "$@"
}

install_cran_pkg() {
  local pkg="$1"
  log "Installing CRAN package in fresh R session: $pkg"

  PKG="$pkg" CRAN_REPO="$CRAN_REPO" Rscript --vanilla - <<'RSCRIPT'
pkg <- Sys.getenv("PKG")
cran_repo <- Sys.getenv("CRAN_REPO")

options(
  repos = c(CRAN = cran_repo),
  timeout = 1000,
  Ncpus = 1
)

cat("\n==== Installing CRAN package:", pkg, "====\n")
cat("CRAN:", cran_repo, "\n")

install.packages(
  pkg,
  dependencies = c("Depends", "Imports", "LinkingTo"),
  Ncpus = 1
)

if (!requireNamespace(pkg, quietly = TRUE)) {
  stop("Package was not loadable after install: ", pkg, call. = FALSE)
}

cat("OK:", pkg, as.character(packageVersion(pkg)), "\n")
RSCRIPT
}

install_local_tarball() {
  local file="$1"
  local path="$PATCH_DIR/$file"

  if [[ ! -f "$path" ]]; then
    log "Local tarball not found, skipping: $path"
    return 0
  fi

  log "Installing local tarball in fresh R session: $file"

  LOCAL_TARBALL="$path" Rscript --vanilla - <<'RSCRIPT'
path <- Sys.getenv("LOCAL_TARBALL")

options(
  repos = c(CRAN = Sys.getenv("CRAN_REPO", "https://packagemanager.posit.co/cran/2022-05-15")),
  timeout = 1000,
  Ncpus = 1
)

cat("\n==== Installing local tarball:", basename(path), "====\n")
install.packages(
  path,
  repos = NULL,
  type = "source",
  dependencies = FALSE,
  Ncpus = 1
)

pkg <- sub("_.*$", "", basename(path))
if (!requireNamespace(pkg, quietly = TRUE)) {
  warning("Package from tarball may not be loadable by guessed name: ", pkg, call. = FALSE)
} else {
  cat("OK:", pkg, as.character(packageVersion(pkg)), "\n")
}
RSCRIPT
}

install_bioc_url() {
  local url="$1"
  local pkg="$2"

  log "Installing Bioconductor package from URL in fresh R session: $pkg"

  BIOC_URL="$url" PKG="$pkg" CRAN_REPO="$CRAN_REPO" Rscript --vanilla - <<'RSCRIPT'
url <- Sys.getenv("BIOC_URL")
pkg <- Sys.getenv("PKG")

options(
  repos = c(CRAN = Sys.getenv("CRAN_REPO")),
  timeout = 1000,
  Ncpus = 1
)

cat("\n==== Installing Bioconductor URL:", pkg, "====\n")
cat(url, "\n")

install.packages(
  url,
  repos = NULL,
  type = "source",
  dependencies = FALSE,
  Ncpus = 1
)

if (!requireNamespace(pkg, quietly = TRUE)) {
  stop("Bioconductor package was not loadable after install: ", pkg, call. = FALSE)
}

cat("OK:", pkg, as.character(packageVersion(pkg)), "\n")
RSCRIPT
}

log "MTD directory: $MTD_DIR"
log "Patch directory: $PATCH_DIR"
log "Log file: $LOG_FILE"
log "Conda env: $ENV_NAME"
log "CONDA_ENSURE: $CONDA_ENSURE"
log "REPAIR_BIOC_CORE: $REPAIR_BIOC_CORE"
log "REPAIR_TIDYVERSE: $REPAIR_TIDYVERSE"
log "INSTALL_OPTIONAL: $INSTALL_OPTIONAL"
log "CRAN_REPO: $CRAN_REPO"

[[ -f "$CONDA_SH" ]] || die "conda.sh not found: $CONDA_SH"

set +u
source "$CONDA_SH"
conda activate "$ENV_NAME"
set -Eeo pipefail

[[ -n "${CONDA_PREFIX:-}" ]] || die "CONDA_PREFIX is empty after conda activate $ENV_NAME"
[[ -x "$CONDA_PREFIX/bin/Rscript" ]] || die "Rscript not found in $CONDA_PREFIX/bin/Rscript"

log "Conda prefix: $CONDA_PREFIX"
log "R: $(command -v R)"
log "Rscript: $(command -v Rscript)"
log "R version: $(R --version | head -n 1)"

R_VERSION="$(Rscript --vanilla -e 'cat(as.character(getRversion()))')"
if [[ "$R_VERSION" != 4.0.* ]]; then
  log "WARNING: this script was tuned for R 4.0.x; detected R $R_VERSION"
fi

cleanup_locks

export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/share/pkgconfig:${PKG_CONFIG_PATH:-}"
export CPATH="$CONDA_PREFIX/include:$CONDA_PREFIX/include/freetype2:${CPATH:-}"
export LIBRARY_PATH="$CONDA_PREFIX/lib:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export CRAN_REPO

if [[ -x "$CONDA_PREFIX/bin/geos-config" ]]; then
  export GEOS_CONFIG="$CONDA_PREFIX/bin/geos-config"
fi

mkdir -p "$HOME/.R"
cat > "$HOME/.R/Makevars" <<'EOF'
CXX11 = x86_64-conda-linux-gnu-c++
CXX14 = x86_64-conda-linux-gnu-c++
CXX17 = x86_64-conda-linux-gnu-c++

CXX11STD = -std=gnu++14
CXX14STD = -std=gnu++14
CXX17STD = -std=gnu++17

CXXFLAGS += -O2
CXX11FLAGS += -O2
CXX14FLAGS += -O2
CXX17FLAGS += -O2
EOF

export R_MAKEVARS_USER="$HOME/.R/Makevars"

if [[ "$CONDA_ENSURE" == "1" ]]; then
  log "Ensuring Conda/system dependencies"
  conda install -y -n "$ENV_NAME" --override-channels -c conda-forge -c bioconda \
    r-biocmanager r-remotes r-rcurl r-curl r-httr \
    r-textshaping r-ragg r-tidyverse r-car r-rstatix r-ggpubr r-plyr \
    r-lattice r-mass r-matrix r-mgcv r-nlme r-survival \
    freetype fontconfig pkg-config harfbuzz fribidi cairo \
    libpng libtiff jpeg libxml2 curl libcurl openssl libuv udunits2 \
    make cmake autoconf automake libtool c-compiler cxx-compiler fortran-compiler
fi

cleanup_locks

# ------------------------------------------------------------
# Phase 1: bootstrap and local rlang/vctrs.
# Important: after this phase, we end the R session.
# ------------------------------------------------------------
run_r "Bootstrap BiocManager/remotes" - <<'RSCRIPT'
options(
  repos = c(CRAN = Sys.getenv("CRAN_REPO")),
  timeout = 1000,
  Ncpus = 1
)

for (p in c("BiocManager", "remotes")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = c("Depends", "Imports", "LinkingTo"), Ncpus = 1)
  }
  cat("OK:", p, as.character(packageVersion(p)), "\n")
}

cat("R:", as.character(getRversion()), "\n")
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
cat("CRAN:", getOption("repos")[["CRAN"]], "\n")
RSCRIPT

# rlang and vctrs are installed first, then we start fresh R processes.
install_local_tarball "rlang_1.1.2.tar.gz"
install_local_tarball "vctrs_0.6.4.tar.gz"

cleanup_locks

# ------------------------------------------------------------
# Phase 2: RCurl and Bioconductor core.
# ------------------------------------------------------------
install_cran_pkg "bitops"
install_cran_pkg "RCurl"

if [[ "$REPAIR_BIOC_CORE" == "1" ]]; then
  install_bioc_url "https://bioconductor.org/packages/3.12/bioc/src/contrib/BiocGenerics_0.36.1.tar.gz" "BiocGenerics"
  install_bioc_url "https://bioconductor.org/packages/3.12/bioc/src/contrib/S4Vectors_0.28.1.tar.gz" "S4Vectors"
  install_bioc_url "https://bioconductor.org/packages/3.12/bioc/src/contrib/IRanges_2.24.1.tar.gz" "IRanges"
  install_bioc_url "https://bioconductor.org/packages/3.12/data/annotation/src/contrib/GenomeInfoDbData_1.2.4.tar.gz" "GenomeInfoDbData"
  install_bioc_url "https://bioconductor.org/packages/3.12/bioc/src/contrib/GenomeInfoDb_1.26.7.tar.gz" "GenomeInfoDb"
fi

cleanup_locks

# ------------------------------------------------------------
# Phase 3: mandatory tidyverse stack.
# Each package is installed in a fresh R session to avoid namespace-cache bugs.
# ------------------------------------------------------------
if [[ "$REPAIR_TIDYVERSE" == "1" ]]; then
  log "Repairing mandatory tidyverse stack from CRAN snapshot $CRAN_REPO"

  # Build dependencies and core tidyverse pieces.
  # Keep rlang/vctrs out of this list because they are intentionally pinned by local tarballs.
  TIDYVERSE_STACK=(
    glue
    cli
    crayon
    lifecycle
    magrittr
    withr
    ellipsis
    R6
    pkgconfig
    utf8
    fansi
    generics
    cpp11
    pillar
    tidyselect
    tibble
    hms
    prettyunits
    progress
    ps
    processx
    callr
    curl
    openssl
    mime
    jsonlite
    httr
    DBI
    blob
    stringi
    stringr
    purrr
    forcats
    dplyr
    tidyr
    ggplot2
    broom
    dbplyr
    dtplyr
    tzdb
    vroom
    readr
    selectr
    xml2
    rvest
    cellranger
    readxl
    haven
    lubridate
    modelr
    rstudioapi
    reprex
    gargle
    googledrive
    googlesheets4
    conflicted
    tidyverse
  )

  for pkg in "${TIDYVERSE_STACK[@]}"; do
    cleanup_locks
    install_cran_pkg "$pkg"
  done
else
  install_cran_pkg "tidyverse"
fi

cleanup_locks

# ------------------------------------------------------------
# Phase 4: original script requirements and local tarballs.
# ------------------------------------------------------------
CORE_CRAN=(
  plyr
  textshaping
  ragg
  car
  rstatix
  ggpubr
)

for pkg in "${CORE_CRAN[@]}"; do
  cleanup_locks
  install_cran_pkg "$pkg"
done

# Remaining local CRAN tarballs from original script.
install_local_tarball "matrixStats_1.3.0.tar.gz"
install_local_tarball "formatR_1.14.tar.gz"
install_local_tarball "lambda.r_1.2.4.tar.gz"
install_local_tarball "futile.options_1.0.1.tar.gz"
install_local_tarball "futile.logger_1.4.3.tar.gz"
install_local_tarball "RColorBrewer_1.1-3.tar.gz"

# UCSC.utils is not native to Bioc 3.12. Keep it optional.
if [[ "$INSTALL_OPTIONAL" == "1" ]]; then
  log "Optional packages requested"
  install_local_tarball "UCSC.utils_1.0.0.tar.gz" || true
  install_cran_pkg "rgeos" || true
fi

cleanup_locks

# ------------------------------------------------------------
# Final validation in a fresh R session.
# ------------------------------------------------------------
run_r "Final fresh R validation" - <<'RSCRIPT'
options(
  repos = c(CRAN = Sys.getenv("CRAN_REPO")),
  timeout = 1000,
  Ncpus = 1
)

pkgs <- c(
  "BiocManager",
  "plyr",
  "rlang",
  "vctrs",
  "BiocGenerics",
  "S4Vectors",
  "IRanges",
  "GenomeInfoDbData",
  "GenomeInfoDb",
  "matrixStats",
  "formatR",
  "lambda.r",
  "futile.options",
  "futile.logger",
  "RColorBrewer",
  "textshaping",
  "ragg",
  "tidyverse",
  "car",
  "rstatix",
  "ggpubr",
  "RCurl"
)

pkg_ok <- function(p) requireNamespace(p, quietly = TRUE)
pkg_ver <- function(p) if (pkg_ok(p)) as.character(packageVersion(p)) else "-"

status <- data.frame(
  Package = pkgs,
  Installed = vapply(pkgs, pkg_ok, logical(1)),
  Version = vapply(pkgs, pkg_ver, character(1)),
  stringsAsFactors = FALSE
)

cat("\n==== Package validation ====\n")
print(status, row.names = FALSE)

cat("\n==== Load test ====\n")
load_status <- vapply(pkgs, function(p) {
  suppressPackageStartupMessages(require(p, character.only = TRUE, quietly = TRUE))
}, logical(1))

print(data.frame(Package = pkgs, Load_OK = load_status), row.names = FALSE)

bad <- unique(c(status$Package[!status$Installed], names(load_status)[!load_status]))
if (length(bad)) {
  stop("Required packages failed final validation: ", paste(bad, collapse = ", "), call. = FALSE)
}

cat("\n==== Success: all required MTD R packages installed/loadable ====\n")
cat("R:", as.character(getRversion()), "\n")
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
RSCRIPT

cleanup_locks

log "Installation finished successfully"
log "Log saved to: $LOG_FILE"

