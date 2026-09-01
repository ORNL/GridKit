# syntax=docker/dockerfile:1.7

ARG BUILDER_IMAGE=gridkit-dev:latest

# Reuse the tested development toolchain to compile GridKit. LLVM and Enzyme
# are build-time dependencies and do not carry into the runtime stage.
FROM ${BUILDER_IMAGE} AS build

USER root

ARG BUILD_JOBS=4

WORKDIR /src
COPY . .

# Enzyme 0.0.206 is known to work with the repository's RelWithDebInfo (-O2)
# configuration. The default Release (-O3) configuration is not used here.
RUN set -eux; \
    cmake -S /src -B /build -G Ninja \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_C_COMPILER=/usr/bin/clang-16 \
      -DCMAKE_CXX_COMPILER=/usr/bin/clang++-16 \
      -DCMAKE_INSTALL_PREFIX=/opt/gridkit \
      -DBUILD_SHARED_LIBS=ON \
      -DBUILD_TESTING=OFF \
      -DGridKit_ENABLE_ENZYME=ON \
      -DGridKit_ENABLE_IPOPT=OFF \
      -DGridKit_ENABLE_SUNDIALS=ON \
      -DENZYME_DIR=/opt/gridkit-deps/enzyme \
      -DSUNDIALS_DIR=/opt/gridkit-deps/sundials \
      -DCMAKE_PREFIX_PATH="/opt/gridkit-deps/sundials;/opt/gridkit-deps/enzyme"; \
    cmake --build /build --parallel "${BUILD_JOBS}"; \
    DESTDIR=/stage cmake --install /build --strip


FROM ubuntu:24.04 AS runtime

RUN set -eux; \
    apt-get -o APT::Sandbox::User=root update; \
    apt-get -o APT::Sandbox::User=root install -y --no-install-recommends \
      libgomp1 \
      libklu2 \
      libstdc++6; \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /stage/opt/gridkit/ /opt/gridkit/
COPY --from=build /opt/gridkit-deps/sundials/lib/ /opt/gridkit/lib/

RUN set -eux; \
    printf '%s\n' /opt/gridkit/lib > /etc/ld.so.conf.d/gridkit.conf; \
    ldconfig; \
    mkdir /work

ENV PATH="/opt/gridkit/bin:${PATH}"

WORKDIR /work

# Keep the image application-neutral. Commands supplied after the image name
# select any installed GridKit executable; an interactive run opens a shell.
CMD ["/bin/bash"]
