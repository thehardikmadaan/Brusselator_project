FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build tools, CMake, and the Eigen3 development library
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libeigen3-dev \
    libvtk9-dev \
    qtbase5-dev \
    libopenmpi-dev \
    libgl1-mesa-dev \
    default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Set working directory inside the container
WORKDIR /app

# Copy the entire project into the container
COPY . .

# Configure, compile, and install the project
RUN mkdir -p build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=.. .. && \
    make -j$(nproc) && \
    make install

# Add the bin directory to the system PATH
ENV PATH="/app/bin:${PATH}"

# Add Java library path for libjawt.so and libjvm.so
ENV LD_LIBRARY_PATH="/usr/lib/jvm/default-java/lib:/usr/lib/jvm/java-11-openjdk-arm64/lib:/usr/lib/jvm/java-11-openjdk-arm64/lib/server:/usr/lib/jvm/java-11-openjdk-amd64/lib:/usr/lib/jvm/java-11-openjdk-amd64/lib/server:${LD_LIBRARY_PATH}"

# Default to bash shell
CMD ["bash"]
