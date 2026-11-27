# Installation

This page explains how to install MTD Explorer from the GitHub repository using the official installation script.

Before starting, review the [system requirements](requirements.md).

!!! warning "Development version"

    MTD Explorer is under active development. Installation commands, database requirements, and output structures may change before the first stable release.

## Functional data resources

Before updating an existing installation, back up the historical functional dictionaries and optional legacy GMT. New installations must supply local GO and KO dictionaries before running the pipeline. See [Local functional resources](functional-resources.md) for acquisition, migration and offline configuration.

## Install Git

Git is required to download the MTD Explorer source code and to keep the installation updated.

Check whether Git is already installed:

```bash
git --version
```

When Git is available, the command displays its installed version.

If Git is not installed, install it using the package manager of your Linux distribution.

### Ubuntu and Linux Mint

```bash
sudo apt update
sudo apt install -y git
```

Confirm the installation:

```bash
git --version
```

## Clone the MTD Explorer repository

After installing Git, download the MTD Explorer repository.

### Using HTTPS

```bash
cd ~
git clone https://github.com/patrick-douglas/MTD-Explorer.git
cd MTD-Explorer
```

HTTPS is the recommended method for users who have not configured an SSH key for GitHub.

### Using SSH

Users with a GitHub SSH key already configured may use:

```bash
cd ~
git clone git@github.com:patrick-douglas/MTD-Explorer.git
cd MTD-Explorer
```

Confirm that the repository was downloaded correctly:

```bash
git status
git log -1 --oneline
```

The installation commands in the following sections should be executed from inside the MTD Explorer repository directory.

## Before running the installer

Before starting the installation, confirm that:

- you are using a supported GNU/Linux system;
- Git is installed;
- the repository was cloned successfully;
- the installation cache is located on a disk with enough free space;
- the selected Miniconda path can be safely created or replaced;
- you have a stable internet connection;
- you have `sudo` access for installing system dependencies;
- the installer is being run from an interactive terminal.

Check the current repository state with:

```bash
git status
git log -1 --oneline
```

Check available disk space with:

```bash
df -h
```

Check memory and CPU resources with:

```bash
free -h
nproc
```

## Installer command

The MTD Explorer installer is executed through `Install.sh`.

Display the current command-line help with:

```bash
bash Install.sh -h
```

The general syntax is:

```text
Install.sh -o <installation-cache> [options]
```

## What the installer does

The installer performs the main steps required to prepare MTD Explorer for use.

In summary, it:

1. validates installer arguments and paths;
2. prepares a persistent installation cache;
3. installs Miniconda when needed;
4. creates the required Conda environments;
5. installs system, Python, R, and bioinformatics dependencies;
6. downloads and validates reference files;
7. prepares HUMAnN databases;
8. prepares Kraken2 taxonomy and microbial databases;
9. prepares the validated Virus-Host DB reference mirror;
10. records installation paths used by MTD Explorer;
11. leaves the installation ready for validation with `MTD_check_installation.sh`.

!!! important "Always verify the installation"

    A completed installer run should always be followed by:

    ```bash
    bash MTD_check_installation.sh --mode full
    ```

    For clean installations, interrupted downloads, or cache transfers between computers, use:

    ```bash
    bash MTD_check_installation.sh --mode deep
    ```

## Required option

### Installation cache

The `-o` option specifies the persistent installation cache:

```text
-o PATH
```

This option is required.

The directory is created and populated automatically when it does not already exist. Files already present in a valid cache may be reused during subsequent or interrupted installations.

Example:

```bash
bash Install.sh \
  -o /home/user/MTD_install_cache
```

The cache may also be located on a separate mounted disk:

```bash
bash Install.sh \
  -o /path/to/large/storage/MTD_install_cache
```

!!! important "Use an absolute path"

    An absolute cache path is recommended.

    Correct:

    ```text
    /home/user/MTD_install_cache
    ```

    Correct:

    ```text
    /media/user/storage/MTD_install_cache
    ```

    Incorrect:

    ```text
    /home/user/~/MTD_install_cache
    ```

    The `~` character is expanded only when it appears at the beginning of a shell path.

## Installer options

| Option | Argument | Description |
|---|---|---|
| `-o` | `PATH` | Persistent installation cache. This option is required. |
| `-p` | `PATH` | Miniconda installation directory. The default is `$HOME/miniconda3`. |
| `-k` | `INT` | Kraken2 k-mer length used when building databases. |
| `-m` | `INT` | Kraken2 minimizer length used when building databases. |
| `-s` | `INT` | Kraken2 minimizer-spaces value used when building databases. |
| `-r` | `INT` | Bracken read length. The default is `75`. |
| `-h` | — | Display the installer help message and exit. |

!!! note "Kraken2 database parameters"

    Most users should leave `-k`, `-m`, and `-s` unset.

    These parameters alter how Kraken2 databases are built and should be changed only when there is a specific technical reason to use non-default database settings.

The installer also accepts the environment variable
`MTD_KRAKEN2_MASKER_THREADS` to control the number of threads used by `k2mask`
during low-complexity masking. When it is not set, the installer uses its
detected CPU thread count (`nproc`).

## Standard installation

For a standard installation using the default Miniconda location and a Bracken read length of 75:

```bash
cd ~/MTD-Explorer

bash Install.sh \
  -o /home/user/MTD_install_cache
```

Replace `/home/user/MTD_install_cache` with the desired cache location.

## Custom Miniconda location

By default, Miniconda is installed at:

```text
$HOME/miniconda3
```

A different location can be selected with `-p`:

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  -p /path/to/miniconda3
```

The selected parent directory must be writable by the current user.

## Existing Miniconda installation

The installer automatically downloads and installs Miniconda.

If the selected Miniconda directory already exists, the installer displays a warning and requests explicit confirmation before permanently deleting that directory.

!!! danger "Existing Miniconda directory"

    Confirming the removal of an existing Miniconda directory permanently deletes the environments and packages stored inside it.

    Before confirming:

    1. verify that the displayed path is correct;
    2. confirm that no important Conda environments are stored there;
    3. back up anything that must be preserved.

    Do not confirm the deletion when the installer displays an unexpected directory.

## Bracken read length

The default Bracken read length is:

```text
75
```

A different value can be selected using `-r`:

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  -r 100
```

The selected value should correspond to the read length intended for Bracken abundance estimation.

Changing this option affects the Bracken files generated during database preparation.

## Save the installation log

The complete installer output should be saved, especially during clean installations or tests on a new computer.

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  2>&1 | tee MTD_installation.log
```

This displays the installation output in the terminal while also saving it to:

```text
MTD_installation.log
```

For a timestamped log:

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  2>&1 | tee "MTD_installation_$(date +%Y%m%d_%H%M%S).log"
```

## During installation

The installer may display status messages such as:

| Status | Meaning |
|---|---|
| `INFO` | Describes the current installation step |
| `PASS` | A required step or validation completed successfully |
| `WARNING` | A non-fatal condition requires attention |
| `RETRY` | A failed operation is being attempted again |
| `ERROR` | A required operation failed |

Large reference databases may require considerable download and processing time.

The installation cache should not be deleted after a successful installation. It can be reused when:

- reinstalling MTD Explorer;
- recovering from an interrupted installation;
- installing the pipeline on another computer;
- validating downloaded reference files;
- avoiding repeated large downloads.

## Dedicated runtime environments

The current installation isolates several tools whose modern dependencies should not replace packages in the legacy main `MTD` environment. In particular, MTD Explorer uses dedicated environments for fastp (`MTD_fastp`), Kraken2/Bracken (`MTD_kraken2`), HUMAnN/MetaPhlAn (`MTD_humann`), and host quantification with featureCounts (`MTD_featurecounts`).

The installation checker currently validates these dedicated production runtimes:

| Environment | Tool/package | Expected version |
|---|---|---:|
| `MTD_fastp` | fastp | `1.3.6` |
| `MTD_featurecounts` | Subread / featureCounts | `2.1.1` |
| `MTD_kraken2` | Kraken2 | `2.17.1` |
| `MTD_kraken2` | Bracken package | `3.1p1` |

The legacy Subread package may remain in the main `MTD` environment for dependency-stack stability, but production host quantification uses the `featureCounts` executable from `MTD_featurecounts`.

Similarly, production Kraken2 and Bracken commands are taken from `MTD_kraken2`. The checker treats mismatches for the featureCounts, Kraken2, and Bracken production runtimes as installation failures; the fastp version check is currently reported as a warning-level version contract while the environment and executable remain required.

!!! note "Bracken package version versus CLI banner"

    The validated Bioconda package is `3.1p1`, but its current upstream
    `bracken -v` banner can still report `v3.0.1`. MTD Explorer therefore
    validates the Conda package metadata for Bracken rather than treating that
    stale executable banner as the package version.

## Kraken2 database cleanup after installation

After the default microbial Kraken2 database and its Bracken distribution are successfully validated, the installer removes the redundant final-installation directory:

```text
kraken2DB_micro/library/bacteria/all
```

The consolidated database library remains available through `library/bacteria/library.fna`. This cleanup reduces the final installation footprint. It does **not** remove the corresponding persistent installation-cache copy, so cached source genomes can still be reused for later installation or rebuild operations.

## Optional HPC node installation

The ordinary installation remains the prerequisite for normal local use. If the pipeline will also run through Slurm, prepare the compute nodes separately after the main installation:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user your_user
```

Then validate the node-local runtime with `Installation/HPC/MTD_check_HPC_nodes.sh`. HPC node installation is optional and is documented separately in [HPC / Slurm execution](../user-guide/hpc-slurm.md).

## After installation

Do not begin a real analysis solely because `Install.sh` reached the end.

First, run the dedicated installation verification procedure described in [Verify the installation](verify-installation.md).
