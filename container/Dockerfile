# syntax=docker/dockerfile:1.7
# Reproduction scaffold for the HPC environment captured on 2026-07-24.
#
# The compute node runs RHEL 8.10. Rocky Linux 8.10 is used as the default
# ABI-compatible public base. Override BASE_IMAGE with an entitled RHEL 8.10
# image when byte-for-byte distribution identity is required.
ARG BASE_IMAGE=rockylinux/rockylinux:8.10@sha256:e8a49c5403b687db05d4d67333fa45808fbe74f36e683cec7abb1f7d0f2338c6
FROM ${BASE_IMAGE}

ARG R_VERSION=4.4.2
ARG RENV_VERSION=1.2.3
ARG PIP_VERSION=25.3
ARG PYTHON_LOCK=requirements-repository-all-target.lock.txt

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN --mount=type=secret,id=build_ca_bundle \
    dnf_ca=() \
    && if [[ -s /run/secrets/build_ca_bundle ]]; then \
        dnf_ca+=(--setopt=sslcacert=/run/secrets/build_ca_bundle); \
    fi \
    && dnf "${dnf_ca[@]}" -y install dnf-plugins-core epel-release \
    && (dnf config-manager --set-enabled powertools || true) \
    && dnf "${dnf_ca[@]}" -y install \
        gcc-toolset-13 \
        gcc-toolset-13-gcc-gfortran \
        make \
        tar \
        gzip \
        bzip2-devel \
        xz-devel \
        zlib-devel \
        pcre2-devel \
        readline-devel \
        libcurl-devel \
        libxml2-devel \
        openssl-devel \
        openblas-devel \
        cairo-devel \
        fontconfig-devel \
        freetype-devel \
        harfbuzz-devel \
        fribidi-devel \
        libpng-devel \
        libjpeg-turbo-devel \
        libtiff-devel \
        libX11-devel \
        libXt-devel \
        tcl-devel \
        tk-devel \
        java-17-openjdk-devel \
        ImageMagick \
        ImageMagick-c++-devel \
        pandoc \
        python39 \
        python39-pip \
        git \
        aria2 \
        curl \
        ca-certificates \
        which \
    && dnf clean all

ENV PATH=/opt/rh/gcc-toolset-13/root/usr/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/rh/gcc-toolset-13/root/usr/lib64

RUN --mount=type=secret,id=build_ca_bundle \
    curl_ca=() \
    && if [[ -s /run/secrets/build_ca_bundle ]]; then \
        curl_ca+=(--cacert /run/secrets/build_ca_bundle); \
    fi \
    && r_downloaded=false \
    && for attempt in {1..8}; do \
        if curl "${curl_ca[@]}" --http1.1 -fsSL \
            "https://cloud.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz" \
            -o "/tmp/R-${R_VERSION}.tar.gz"; then \
            r_downloaded=true; \
            break; \
        fi; \
        sleep 2; \
    done \
    && [[ "${r_downloaded}" == true ]] \
    && tar -C /tmp -xzf "/tmp/R-${R_VERSION}.tar.gz" \
    && cd "/tmp/R-${R_VERSION}" \
    && ./configure \
        --prefix=/opt/R/${R_VERSION} \
        --enable-R-shlib \
        --with-x \
        --with-blas \
        --with-lapack \
    && make -j"$(nproc)" \
    && make install \
    && rm -rf "/tmp/R-${R_VERSION}" "/tmp/R-${R_VERSION}.tar.gz"

ENV PATH=/opt/R/${R_VERSION}/bin:${PATH}
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk
ENV RENV_VERSION=${RENV_VERSION}

# The HPC R module exposes cURL/libcurl 8.7.1. Rocky Linux 8 ships libcurl
# 7.61.1, which is too old to compile the captured R curl 7.1.0 package.
# Keep the distro libcurl for RPM compatibility and install the HPC-matching
# cURL in an isolated prefix.
RUN --mount=type=secret,id=build_ca_bundle \
    dnf_ca=() \
    && if [[ -s /run/secrets/build_ca_bundle ]]; then \
        dnf_ca+=(--setopt=sslcacert=/run/secrets/build_ca_bundle); \
    fi \
    && dnf "${dnf_ca[@]}" -y install cmake \
    && dnf clean all

ARG CURL_VERSION=8.7.1
ARG CURL_SHA256=6fea2aac6a4610fbd0400afb0bcddbe7258a64c63f1f68e5855ebc0c659710cd

RUN --mount=type=secret,id=build_ca_bundle \
    curl_ca=() \
    && if [[ -s /run/secrets/build_ca_bundle ]]; then \
        curl_ca+=(--cacert /run/secrets/build_ca_bundle); \
    fi \
    && curl_downloaded=false \
    && for attempt in {1..8}; do \
        if curl "${curl_ca[@]}" --http1.1 -fsSL \
            "https://curl.se/download/curl-${CURL_VERSION}.tar.xz" \
            -o "/tmp/curl-${CURL_VERSION}.tar.xz"; then \
            curl_downloaded=true; \
            break; \
        fi; \
        sleep 2; \
    done \
    && [[ "${curl_downloaded}" == true ]] \
    && printf '%s  %s\n' "${CURL_SHA256}" "/tmp/curl-${CURL_VERSION}.tar.xz" \
        | sha256sum --check --strict \
    && tar -C /tmp -xJf "/tmp/curl-${CURL_VERSION}.tar.xz" \
    && cd "/tmp/curl-${CURL_VERSION}" \
    && ./configure \
        --prefix="/opt/curl/${CURL_VERSION}" \
        --with-openssl \
        --with-zlib \
        --with-ca-bundle=/etc/pki/tls/certs/ca-bundle.crt \
        --disable-static \
        --enable-shared \
    && make -j"$(nproc)" \
    && make install \
    && printf '%s\n' "/opt/curl/${CURL_VERSION}/lib" \
        > "/etc/ld.so.conf.d/curl-${CURL_VERSION}.conf" \
    && ldconfig \
    && rm -rf "/tmp/curl-${CURL_VERSION}" "/tmp/curl-${CURL_VERSION}.tar.xz"

ENV CURL_VERSION=${CURL_VERSION}
ENV PATH=/opt/curl/${CURL_VERSION}/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/curl/${CURL_VERSION}/lib:${LD_LIBRARY_PATH}
ENV PKG_CONFIG_PATH=/opt/curl/${CURL_VERSION}/lib/pkgconfig

WORKDIR /opt/soft-coupling-environment
COPY locks/ /opt/soft-coupling-environment/locks/
COPY scripts/verify_environment.R /opt/soft-coupling-environment/scripts/verify_environment.R
COPY scripts/verify_python_environment.py /opt/soft-coupling-environment/scripts/verify_python_environment.py

RUN --mount=type=secret,id=build_ca_bundle \
    if [[ -s /run/secrets/build_ca_bundle ]]; then \
        export CURL_CA_BUNDLE=/run/secrets/build_ca_bundle \
            SSL_CERT_FILE=/run/secrets/build_ca_bundle; \
    fi \
    && Rscript -e "install.packages(sprintf('https://cloud.r-project.org/src/contrib/renv_%s.tar.gz', Sys.getenv('RENV_VERSION')), repos = NULL, type = 'source')" \
    && Rscript -e "renv::consent(provided = TRUE); renv::activate(project = '.'); lock <- renv::lockfile_read('locks/renv.lock'); packages <- names(lock\$Packages)[vapply(lock\$Packages, function(record) !identical(record\$Source, 'R'), logical(1L))]; renv::restore(lockfile = 'locks/renv.lock', packages = packages, prompt = FALSE)"

# The captured HPC library has a newer recommended cluster package than the
# version bundled with the R 4.4.2 source tarball. It is not pulled transitively
# by renv, so install that exact target version into the R installation library.
# Disable the activation profile written by renv: restore installed these locked
# packages into R's installation library, while that profile points a later R
# process at an empty project library and hides them.
ENV R_PROFILE_USER=/dev/null

RUN --mount=type=secret,id=build_ca_bundle \
    curl_ca=() \
    && if [[ -s /run/secrets/build_ca_bundle ]]; then \
        curl_ca+=(--cacert /run/secrets/build_ca_bundle); \
    fi \
    && cluster_downloaded=false \
    && for attempt in {1..8}; do \
        if curl "${curl_ca[@]}" --http1.1 -fsSL \
            "https://cloud.r-project.org/src/contrib/Archive/cluster/cluster_2.1.8.1.tar.gz" \
            -o /tmp/cluster_2.1.8.1.tar.gz; then \
            cluster_downloaded=true; \
            break; \
        fi; \
        sleep 2; \
    done \
    && [[ "${cluster_downloaded}" == true ]] \
    && R CMD INSTALL /tmp/cluster_2.1.8.1.tar.gz \
    && rm -f /tmp/cluster_2.1.8.1.tar.gz \
    && Rscript --vanilla -e "stopifnot(identical(as.character(packageVersion('cluster')), '2.1.8.1'))"

RUN --mount=type=secret,id=build_ca_bundle \
    build_pip_ca=/tmp/soft-coupling-build-ca-bundle.pem \
    && trap 'rm -f /tmp/soft-coupling-build-ca-bundle.pem' EXIT \
    && if [[ -s /run/secrets/build_ca_bundle ]]; then \
        cp /etc/pki/tls/certs/ca-bundle.crt "${build_pip_ca}"; \
        cat /run/secrets/build_ca_bundle >> "${build_pip_ca}"; \
        chmod 0600 "${build_pip_ca}"; \
        export PIP_CERT="${build_pip_ca}" \
            REQUESTS_CA_BUNDLE="${build_pip_ca}"; \
    fi \
    && python3.9 -m pip install --upgrade "pip==${PIP_VERSION}" \
    && python3.9 -m pip install --require-hashes -r "locks/${PYTHON_LOCK}" \
    && Rscript scripts/verify_environment.R locks/packages.lock.tsv target \
    && python3.9 scripts/verify_python_environment.py "locks/${PYTHON_LOCK}" \
    && rm -f "${build_pip_ca}" \
    && test ! -e "${build_pip_ca}"

ARG IMAGEMAGICK7_VERSION=7.1.2-29
ARG IMAGEMAGICK7_SHA256=4b131411aa77b051908ea35df65ca2c52dbe04b48c3ef6b06f35c1bb595b212f

# Keep Rocky's ImageMagick 6 libraries for the already-verified R magick
# package, and add an isolated, digest-verified ImageMagick 7 installation for
# the unified `magick` CLI. The limited policy is safer for a public image than
# the upstream open default.
RUN --mount=type=secret,id=build_ca_bundle \
    curl_ca=() \
    && if [[ -s /run/secrets/build_ca_bundle ]]; then \
        curl_ca+=(--cacert /run/secrets/build_ca_bundle); \
    fi \
    && imagemagick7_downloaded=false \
    && for attempt in {1..8}; do \
        if curl "${curl_ca[@]}" --http1.1 -fsSL \
            "https://github.com/ImageMagick/ImageMagick/releases/download/${IMAGEMAGICK7_VERSION}/ImageMagick-${IMAGEMAGICK7_VERSION}.tar.xz" \
            -o "/tmp/ImageMagick-${IMAGEMAGICK7_VERSION}.tar.xz"; then \
            imagemagick7_downloaded=true; \
            break; \
        fi; \
        sleep 2; \
    done \
    && [[ "${imagemagick7_downloaded}" == true ]] \
    && printf '%s  %s\n' \
        "${IMAGEMAGICK7_SHA256}" \
        "/tmp/ImageMagick-${IMAGEMAGICK7_VERSION}.tar.xz" \
        | sha256sum --check --strict \
    && tar -C /tmp -xJf "/tmp/ImageMagick-${IMAGEMAGICK7_VERSION}.tar.xz" \
    && cd "/tmp/ImageMagick-${IMAGEMAGICK7_VERSION}" \
    && ./configure \
        --prefix="/opt/ImageMagick/${IMAGEMAGICK7_VERSION}" \
        --with-quantum-depth=16 \
        --enable-shared \
        --disable-static \
        --without-perl \
        --with-security-policy=limited \
    && make -j"$(nproc)" \
    && make install \
    && printf '%s\n' "/opt/ImageMagick/${IMAGEMAGICK7_VERSION}/lib" \
        > /etc/ld.so.conf.d/imagemagick7.conf \
    && ldconfig \
    && rm -rf \
        "/tmp/ImageMagick-${IMAGEMAGICK7_VERSION}" \
        "/tmp/ImageMagick-${IMAGEMAGICK7_VERSION}.tar.xz"

ENV IMAGEMAGICK7_VERSION=${IMAGEMAGICK7_VERSION}
ENV PATH=/opt/ImageMagick/${IMAGEMAGICK7_VERSION}/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/ImageMagick/${IMAGEMAGICK7_VERSION}/lib:${LD_LIBRARY_PATH}
ENV PKG_CONFIG_PATH=/opt/ImageMagick/${IMAGEMAGICK7_VERSION}/lib/pkgconfig:${PKG_CONFIG_PATH}

RUN magick -version \
    && magick identify -version \
    && magick -list policy \
    && magick -size 1x1 xc:white /tmp/imagemagick7-smoke.png \
    && magick identify /tmp/imagemagick7-smoke.png \
    && rm -f /tmp/imagemagick7-smoke.png \
    && Rscript --vanilla -e "stopifnot(requireNamespace('magick', quietly = TRUE)); print(magick::magick_config())" \
    && Rscript --vanilla -e "stopifnot(identical(as.character(packageVersion('shadowtext')), '0.1.6')); stopifnot(is.function(shadowtext::shadowtextGrob)); stopifnot(identical(as.character(packageVersion('textshaping')), '1.0.5'))" \
    && Rscript --vanilla -e "stopifnot(identical(as.character(packageVersion('svglite')), '2.2.2')); output <- tempfile(fileext = '.svg'); svglite::svglite(output); plot.new(); text(0.5, 0.5, 'svglite smoke test'); invisible(dev.off()); stopifnot(file.exists(output), file.info(output)\$size > 0); unlink(output)"

WORKDIR /workspace
