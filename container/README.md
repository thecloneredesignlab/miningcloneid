# soft_coupling HPC runtime inventory

This directory records the environment needed by the current
`oxygen/code/O2_supply_demand_MAP` model, all code under its `analysis/`
directory, the transitive O2 simulation/visualization/report/runner code, and
all Python scripts in this repository.

The accepted snapshot was collected by Slurm job `19623212` on compute node
`hpctpa3pc0007.cm.cluster` at `2026-07-24T17:13:25Z`. It corresponds to local
branch `soft_coupling` at commit
`fc4aef4ee32f1c72cdfe55b2f1c74b41aa07dd0d`, including the working-tree code
present during the audit.

## Captured platform

| Component | Captured value |
| --- | --- |
| Operating system | Red Hat Enterprise Linux 8.10 (Ootpa) |
| Architecture | x86_64 |
| Kernel | 4.18.0-553.70.1.el8_10.x86_64 |
| glibc | 2.28 |
| Requested module | `ml R/4.4` |
| Resolved module | `R/4.4.2-gfbf-2024a` |
| R | 4.4.2, `x86_64-pc-linux-gnu` |
| Compiler stack | GCC/G++/gfortran 13.3.0 |
| BLAS stack | OpenBLAS 0.3.27 + FlexiBLAS 3.4.4 |
| Python | `/usr/bin/python3`, CPython 3.9.20 |
| Java | OpenJDK 17.0.6 |

The R module loads 57 modules. The complete ordered stack is in
`manifests/hpc-loaded-modules.tsv`, and the raw module output is retained under
`hpc_snapshot/modules/`.

R used these libraries in priority order:

1. `/home/4482173/R/x86_64-pc-linux-gnu-library/4.4`
2. `/app/eb/software/R/4.4.2-gfbf-2024a/lib64/R/library`

The snapshot contains 759 installed package records across both libraries and
693 unique selected package names. Duplicate names from lower-priority
libraries are retained rather than discarded.

## Code and dependency coverage

The static audit covers:

- 127 R files and one Python wrapper under
  `oxygen/code/O2_supply_demand_MAP/analysis/`.
- 341 R/Rmd/qmd files across the audited O2 and test surfaces.
- 20 Python files across the repository.
- R namespace calls, package loaders, source calls, Rscript dispatch, and
  external-program calls.

`packages.lock.tsv` contains 106 exact R target records: the 95-package O2
runtime plus explicit Docker targets (`magick`, `shadowtext`, `svglite`, and
`textshaping`) and their complete runtime dependency closure. `shadowtext`
retains the current CRAN/HPC-observed `0.1.6` version; `svglite` is fixed to
the current CRAN `2.2.2`; `textshaping` is fixed to the current
CRAN/HPC-observed `1.0.5`; and `systemfonts` retains the HPC-observed
`1.3.2` version. The analysis-only recursive runtime closure contains 40
unique R packages. Base and recommended packages are included.

The exact 95-package HPC-observed lock was independently checked on the same
module stack by Slurm job `19623467`; every installed DESCRIPTION version
matched.

Important direct analysis versions include:

| Package | Version |
| --- | --- |
| cluster | 2.1.8.1 |
| data.table | 1.18.4 |
| dplyr | 1.2.1 |
| Rtsne | 0.17 |
| uwot | 0.2.4 |
| Matrix | 1.7-3 |
| Rcpp | 1.1.1-1.1 |
| shadowtext | 0.1.6 |
| svglite | 2.2.2 |
| textshaping | 1.0.5 |
| systemfonts | 1.3.2 |

The full per-file mapping is
`manifests/r-package-usage-by-file.tsv`. Recursive runtime and optional
dependency closures are in `hpc_snapshot/r/`.

## Observed HPC gaps

The capture deliberately distinguishes what is installed on HPC from what the
repository code requires:

- R package `magick` is referenced by current report code but is not installed
  in the captured HPC libraries. The Docker target pins CRAN `magick` 2.9.1.
- The default HPC Python has only `pip` and `setuptools`; `openpyxl` is missing.
- Repository-wide Python imports also require NumPy, SciPy, Matplotlib and
  PyYAML, none of which are installed in the default HPC Python.
- ImageMagick's `magick` executable and Pandoc are absent from the captured
  HPC PATH. The Docker target adds an isolated ImageMagick 7 installation for
  the `magick` CLI while retaining the system ImageMagick 6 libraries used by
  the locked R `magick` package.
- Five unresolved shared libraries were found in unrelated installed R
  packages (`units`, `BPCells`, and `sf`). None of those packages belongs to
  the O2 or analysis runtime closure; all shared libraries used by the locked
  O2 package set resolved successfully.

See `manifests/environment-gaps.tsv` for the machine-readable gap list.

## Python locks

- `locks/requirements-hpc-observed.txt` records the actual default HPC Python.
- `locks/requirements-o2-target.lock.txt` contains the O2 Python environment:
  `openpyxl==3.1.5` and `et_xmlfile==2.0.0`.
- `locks/requirements-repository-all-target.lock.txt` contains all 17 resolved
  direct and transitive distributions needed by repository Python scripts.

The target locks were resolved for CPython 3.9, Linux x86_64,
manylinux2014, and include wheel SHA-256 hashes. The pip JSON resolution reports
are preserved beside them. Because pip evaluates dependency environment markers
against the resolver interpreter even when `--python-version 3.9` selects target
wheels, the Python 3.9-only `importlib-resources` and `zipp` dependencies are
explicit resolution inputs. The lock generator rejects a repository-wide report
that omits either package.

## Directory guide

- `hpc_snapshot/`: raw compute-node facts, installed packages, package
  dependency edges, system RPMs, shared-library linkage, module outputs, and
  Slurm provenance.
- `manifests/`: code-to-package mappings, code hashes, runtime library mapping,
  analysis/Python dependency closures, environment gaps, and concise summaries.
- `locks/`: compact R and Python installation locks.
- `scripts/`: repeatable capture, lock generation, and verification programs.
- `Dockerfile`: the verified RHEL-compatible reconstruction recipe.

## Docker status

The HPC system is RHEL 8.10, while the default public base is the digest-pinned
`rockylinux/rockylinux:8.10` image. Rocky Linux is ABI-compatible but not
distribution-identical. Use an entitled RHEL 8.10 image through
`--build-arg BASE_IMAGE=...` if exact distribution identity is required.

The HPC R environment is an EasyBuild stack with 57 modules. The Dockerfile
reconstructs the runtime using OS packages and GCC Toolset 13; it does not
claim that those binaries are byte-identical to EasyBuild.

The last published full repository-Python `r44` artifact was built and
verified locally on 2026-07-30. Its recorded digests predate the
`svglite`/`textshaping` target refresh above; rebuild and republish before
treating those image digests as verification of the current R lock:

| Field | Verified value |
| --- | --- |
| Tag | `o2_supply_demand_map:r44` |
| OCI index digest | `sha256:32c49db0ad27a0b5832b601ba96e2b72bfc1e2f1ccbf34687f8f596f1f7cdcd5` |
| linux/amd64 manifest digest | `sha256:6c75181a01a41436251eb63250efb80864b7dee38196bda2729446fe86423327` |
| Platform | `linux/amd64` |
| Size | 1,420,759,130 bytes |
| Container OS | Rocky Linux 8.10 |
| R | 4.4.2 |
| Python | 3.9.25 |
| pip | 25.3 |
| Git | 2.43.7 |
| aria2 | 1.35.0 |
| cURL/libcurl | 8.7.1 |
| ImageMagick CLI | 7.1.2-29, Q16-HDRI, limited security policy |
| R `magick` backend | ImageMagick 6.9.13.50 |
| Locked R records | 106 checked against the pre-refresh lock |
| Locked Python distributions | 17 checked, 0 mismatches |

Build the full repository-Python target with:

```bash
docker buildx build --platform linux/amd64 --load \
  -t o2_supply_demand_map:r44 \
  -f container/Dockerfile container
```

If the build network uses TLS interception, pass a CA-only bundle from outside
the repository as an ephemeral BuildKit secret:

```bash
docker buildx build --platform linux/amd64 --load \
  --secret id=build_ca_bundle,src=/absolute/path/outside/repository/ca-bundle.pem \
  -t o2_supply_demand_map:r44 \
  -f container/Dockerfile container
```

Do not place that bundle in `container/` and do not pass private keys or client
certificates. The secret is mounted only for network operations. The temporary
pip CA merge is explicitly deleted in the same build layer and its absence is
asserted before the layer succeeds. The verified image contains neither the
build secret nor the enterprise CA certificates used for the local build.

Use the smaller O2-only Python lock with:

```bash
docker buildx build --platform linux/amd64 --load \
  --build-arg PYTHON_LOCK=requirements-o2-target.lock.txt \
  -t o2_supply_demand_map:r44-o2 \
  -f container/Dockerfile container
```

## Re-run the audit

Generate static local manifests:

```bash
Rscript container/scripts/audit_code_dependencies.R "$PWD" container/manifests
python3 container/scripts/audit_python_dependencies.py "$PWD" container/manifests
```

Copy `container/scripts/`, `container/manifests/`, and the Python source
snapshot to a shared HPC audit directory, then submit
`capture_hpc_environment.sh`. The accepted capture uses the exact module
initialization implemented by the production submitter:

```bash
source /etc/profile.d/modules.sh
module use /app/eb/modules/all
ml R/4.4
```

After copying the accepted `hpc_snapshot/` back locally, regenerate compact
locks and checksums:

```bash
python3 container/scripts/build_lockfiles.py container
```

## Verification

Verify the captured HPC package versions:

```bash
Rscript container/scripts/verify_environment.R \
  container/locks/packages.lock.tsv hpc-observed
```

Verify a completed Docker target:

```bash
docker run --rm --platform linux/amd64 o2_supply_demand_map:r44 \
  Rscript /opt/soft-coupling-environment/scripts/verify_environment.R \
    /opt/soft-coupling-environment/locks/packages.lock.tsv target

docker run --rm --platform linux/amd64 o2_supply_demand_map:r44 \
  python3.9 /opt/soft-coupling-environment/scripts/verify_python_environment.py \
    /opt/soft-coupling-environment/locks/requirements-repository-all-target.lock.txt

docker run --rm --platform linux/amd64 o2_supply_demand_map:r44 \
  bash -lc 'magick -version &&
    magick -size 2x2 xc:white /tmp/imagemagick7-smoke.png &&
    magick identify /tmp/imagemagick7-smoke.png'
```

`SHA256SUMS` covers all files in this directory except itself. The raw HPC
snapshot has its own independently generated `hpc_snapshot/SHA256SUMS`.
