FROM docker.io/library/julia:1.10.10-bookworm

# Set the working directory.
WORKDIR /app
ENV JULIA_DEPOT_PATH=/app/.julia

# Copy the project files to the container.
COPY . /app

RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()' \
    && chown -R 65534:65534 /app

USER 65534

CMD ["julia", "-e", "using Pkg; Pkg.activate(\".\"); using OpenSEFA; OpenSEFA.start_server()"]

