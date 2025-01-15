FROM compss/compss-tutorial:3.3
LABEL org.opencontainers.image.authors="support-compss@bsc.es"

COPY . Qdislib/

ENV PYTHONPATH=$PYTHONPATH:/Qdislib:/opt/COMPSs/Bindings/python/3/
ENV LC_ALL=C.UTF-8
RUN python3 -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org --upgrade /Qdislib/

ENV COMPSS_LOAD_SOURCE false

# Expose SSH port and run SSHD
EXPOSE 22
CMD ["/usr/sbin/sshd","-D"]
