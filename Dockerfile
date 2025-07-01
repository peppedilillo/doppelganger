FROM lsstsqre/centos:7-stack-lsst_distrib-d_latest

# Create a wrapper script for Python
RUN echo '#!/bin/bash' > /opt/lsst/software/stack/python-lsst && \
    echo 'source /opt/lsst/software/stack/loadLSST.bash' >> /opt/lsst/software/stack/python-lsst && \
    echo 'setup lsst_distrib' >> /opt/lsst/software/stack/python-lsst && \
    echo 'exec python "$@"' >> /opt/lsst/software/stack/python-lsst && \
    chmod +x /opt/lsst/software/stack/python-lsst
