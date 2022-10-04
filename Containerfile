FROM docker.io/ensemblorg/ensembl-vep:release_104.3

# Root needed to install packages; the entrypoint switches to a regular user
USER root

# Install all dependencies
COPY ./scripts/install_dependencies.sh /opt/
RUN bash /opt/install_dependencies.sh

# Install VEP plugins
# --NO_UPDATE since VEP will otherwise abort when a new VEP release is available
RUN ./INSTALL.pl --AUTO p --PLUGINS all --CACHEDIR /opt/vep-plugins/ --NO_UPDATE

# Install VEP plugin LOFTEE (https://github.com/konradjk/loftee)
COPY ./binaries/loftee-v1.0.3.tar.gz /opt/vep-plugins
RUN cd /opt/vep-plugins && \
    tar xvzf ./loftee-v1.0.3.tar.gz && \
    cp -a loftee-1.0.3/* Plugins/

RUN apt-get update && apt-get install python3-ruamel.yaml
RUN pip3 install aush==0.1.3 isal==0.11.1 --no-cache

# Workaround for leak in built-in version
RUN apt-get remove -y python3-pysam
RUN pip3 install pysam==0.16.0.1 --no-cache

# Required for merging VCFs during setup
RUN apt-get install -y vcftools

# Create folder for mounting the (shared) cache
RUN mkdir -p /data/cache && touch /data/cache/not_mounted
# Create folder for user data (i.e. the user's current working directory)
RUN mkdir -p /data/user && touch /data/user/not_mounted

COPY ./pipeline/ /opt/annovep/pipeline/
COPY ./scripts/ /opt/annovep/scripts/
COPY ./annotations/ /opt/annovep/annotations/

# Normalize permissions
RUN find /opt/annovep/ -type f -exec chmod +r \{\} \; && \
    find /opt/annovep/ -type d -exec chmod +rx \{\} \;

# Mountpoint for the current working directory
WORKDIR /data/user

ENTRYPOINT [ "/usr/bin/python3", "/opt/annovep/scripts/entrypoint.py" ]
