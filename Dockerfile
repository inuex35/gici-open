FROM ubuntu:20.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Set up timezone
RUN apt-get update && apt-get install -y tzdata
RUN ln -fs /usr/share/zoneinfo/Asia/Tokyo /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata

# Install build dependencies and required libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libopencv-dev \
    libeigen3-dev \
    libyaml-cpp-dev \
    libgoogle-glog-dev \
    libgflags-dev \
    libatlas-base-dev \
    libsuitesparse-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Ceres Solver
RUN apt-get update && apt-get install -y \
    libgoogle-glog-dev \
    libgflags-dev \
    libatlas-base-dev \
    libsuitesparse-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
    
RUN git clone https://github.com/ceres-solver/ceres-solver.git /opt/ceres-solver && \
    cd /opt/ceres-solver && \
    git checkout 2.1.0 && \
    mkdir build && cd build && \
    cmake .. -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF && \
    make -j$(nproc) && \
    make install

# Create app directory
WORKDIR /app

# Clone GICI-LIB repository
RUN git clone https://github.com/inuex35/gici-open.git . && \
    git checkout forppc2024 && \
    git submodule update --init --recursive

# Build GICI-LIB
RUN mkdir -p build && \
    cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc)

# Set the entrypoint
ENTRYPOINT ["/app/build/gici_main"]
CMD ["./option/tc.yaml"]
