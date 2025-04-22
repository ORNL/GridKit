FROM ubuntu:24.04

RUN apt update -y && apt install -y \
	gcc \
	g++ \
	cmake \
	ninja-build \
	git \
	libopenblas-dev \
	&& apt clean && rm -rf /var/lib/apt/lists/*

# Install SuiteSparse for klu
RUN git clone --single-branch --depth 1 --branch v7.10.2 https://github.com/DrTimothyAldenDavis/SuiteSparse.git /tmp/suitesparse && \
	cd /tmp/suitesparse && \
	mkdir -p build && cd build && \
	cmake -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DSUITESPARSE_ENABLE_PROJECTS="klu" .. && \
	ninja -j $(nproc) && \
	ninja install && \
	rm -rf /tmp/suitesparse

# Install sundials using klu from previous step
RUN git clone --single-branch --depth 1 --branch v7.3.0 https://github.com/LLNL/sundials.git /tmp/sundials && \
	cd /tmp/sundials && \
	mkdir -p build && cd build && \
	cmake -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DENABLE_KLU=ON \
		-DKLU_ROOT=/usr/local/lib/cmake/KLU .. && \
	ninja -j $(nproc) && \
	ninja install && \
	rm -rf /tmp/sundials
