FROM docker.io/ensemblorg/ensembl-vep:release_104.3

# Root needed to install packages; the entrypoint switches to a regular user
USER root

# Install VEP plugins
# --NO_UPDATE since VEP will otherwise abort when a new VEP release is available
RUN ./INSTALL.pl --AUTO p --PLUGINS all --CACHEDIR /opt/vep-plugins/ --NO_UPDATE

# Install VEP plugin LOFTEE (https://github.com/konradjk/loftee)
COPY ./binaries/loftee-v1.0.3.tar.gz /opt/vep-plugins
RUN cd /opt/vep-plugins && \
    tar xvzf ./loftee-v1.0.3.tar.gz && \
    cp -a loftee-1.0.3/* Plugins/

# Install dependencies for pipeline scripts
COPY ./scripts/install_dependencies.sh /opt/
RUN bash /opt/install_dependencies.sh

# Pip must be upgraded in order to enable installation of requirements
RUN python3.7 -m pip install --upgrade pip==23.3.1

# Install specific python dependencies for containerized version
COPY requirements.container.txt /opt/annovep/
RUN python3.7 -m pip install --no-cache -r /opt/annovep/requirements.container.txt

# Create folder for mounting the (shared) cache
RUN mkdir -p /data/cache && touch /data/cache/not_mounted
# Create folder for user data (i.e. the user's current working directory)
RUN mkdir -p /data/user && touch /data/user/not_mounted

COPY ./annovep /opt/annovep/annovep/
COPY ./MANIFEST.in /opt/annovep/
COPY ./scripts/ /opt/annovep/scripts/
COPY ./setup.py /opt/annovep/
RUN python3.7 -m pip install /opt/annovep/

# Normalize permissions
RUN find /opt/annovep/ -type f -exec chmod +r \{\} \; && \
    find /opt/annovep/ -type d -exec chmod +rx \{\} \;

# Mountpoint for the current working directory
WORKDIR /data/user

ENTRYPOINT [ "/usr/bin/python3.7", "/opt/annovep/scripts/entrypoint.py" ]
