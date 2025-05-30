FROM ubuntu:24.04

RUN apt update -y && apt install -y --no-install-recommends\
	gcc \
	g++ \
	make \
	cmake \
	ninja-build \
	git \
	ca-certificates \
	libopenblas-dev \
	&& apt clean && rm -rf /var/lib/apt/lists/*

# Install SuiteSparse for klu
RUN git clone --single-branch --depth 1 --branch v7.10.2 https://github.com/DrTimothyAldenDavis/SuiteSparse.git /tmp/suitesparse && \
	cd /tmp/suitesparse && \
	mkdir -p build && cd build && \
	cmake -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DSUITESPARSE_ENABLE_PROJECTS="klu" .. && \
	ninja install && \
	rm -rf /tmp/suitesparse

# Install sundials using klu from previous step
RUN git clone --single-branch --branch develop https://github.com/LLNL/sundials.git /tmp/sundials && \
	cd /tmp/sundials && \
	git checkout 07d21c2c5ae33211a9a2fafd8eac56de5582dce0 && \
	mkdir -p build && cd build && \
	cmake -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DENABLE_KLU=ON .. && \
	ninja install && \
	rm -rf /tmp/sundials
