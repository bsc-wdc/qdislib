# Binder image for the Qdislib notebooks (mybinder.org uses binder/Dockerfile automatically).
#
# Base: the official COMPSs tutorial image, which already ships PyCOMPSs, Qiskit, Qiskit Aer,
# Qibo (>=0.3.5, needed for NumPy 2), matplotlib and JupyterLab. Pinned by digest (multi-arch
# index, 2026-09-25) so the demo cannot break if the "latest" tag changes; update it on purpose.
FROM compss/compss-tutorial:latest@sha256:7467e9f320d5d1de7f073f9590d2ce429b634034b710e72c3767bdb4c555b213

# Binder runs the container as UID 1000; in this image that user is "ubuntu". mybinder.org
# builds with "--build-arg NB_USER=jovyan", a user this image does not have, so the name is
# fixed here and NB_USER is accepted but not used.
ARG NB_USER
ARG NB_UID=1000
ENV USER=ubuntu \
    NB_UID=${NB_UID} \
    HOME=/home/ubuntu

USER root

# COMPSs starts its local worker with "ssh -p 22 localhost ...". Binder users cannot run a
# server on port 22, so a private sshd listens on 127.0.0.1:2222 (see start.sh) and this
# wrapper sends COMPSs' localhost connections there (COMPSs passes "-p 22" as one argument).
# Other ssh calls are left untouched.
RUN printf '%s\n' '#!/bin/bash' \
      'args=("$@")' \
      'if [[ " $* " == *" localhost "* || " $* " == *" 127.0.0.1 "* ]]; then' \
      '  for i in "${!args[@]}"; do' \
      '    if [[ "${args[$i]}" =~ ^-p[[:space:]]*22$ ]]; then args[$i]="-p2222"; fi' \
      '    if [[ "${args[$i]}" == "-p" && "${args[$((i+1))]}" == "22" ]]; then args[$((i+1))]=2222; fi' \
      '  done' \
      'fi' \
      'exec /usr/bin/ssh "${args[@]}"' > /usr/local/bin/ssh \
 && chmod 755 /usr/local/bin/ssh

COPY . ${HOME}
# Install Qdislib from this repository, so the notebooks always run the code they ship with.
# graphviz (Python package) draws the PyCOMPSs task graph; the "dot" binary is in the base image.
# [cache] and ipywidgets are for the other notebooks; Qiskit and Qibo are in the base image.
RUN python3 -m pip install --no-cache-dir --break-system-packages "${HOME}[cache]" graphviz ipywidgets \
 && chown -R ${NB_UID}:${NB_UID} ${HOME}

USER ubuntu
WORKDIR ${HOME}

# Passwordless SSH to the private sshd. PermitUserEnvironment passes the COMPSs environment
# (taken from a login shell, i.e. /etc/profile.d/compss.sh) into the worker's SSH session.
RUN mkdir -p ${HOME}/.ssh ${HOME}/.sshd && chmod 700 ${HOME}/.ssh \
 && ssh-keygen -q -t ed25519 -N "" -f ${HOME}/.ssh/id_ed25519 \
 && cp ${HOME}/.ssh/id_ed25519.pub ${HOME}/.ssh/authorized_keys \
 && ssh-keygen -q -t ed25519 -N "" -f ${HOME}/.sshd/ssh_host_ed25519_key \
 && printf '%s\n' "Port 2222" "ListenAddress 127.0.0.1" \
      "HostKey ${HOME}/.sshd/ssh_host_ed25519_key" "PidFile ${HOME}/.sshd/sshd.pid" \
      "AuthorizedKeysFile .ssh/authorized_keys" "PasswordAuthentication no" "UsePAM no" \
      "StrictModes no" "PermitUserEnvironment yes" > ${HOME}/.sshd/sshd_config \
 && printf '%s\n' "Host localhost 127.0.0.1" "  StrictHostKeyChecking no" \
      "  UserKnownHostsFile /dev/null" "  LogLevel ERROR" > ${HOME}/.ssh/config \
 && bash -lc env | grep -E '^(PATH|CLASSPATH|LD_LIBRARY_PATH|DLB_HOME|JAVA_HOME)=' > ${HOME}/.ssh/environment \
 && chmod 600 ${HOME}/.ssh/config ${HOME}/.ssh/authorized_keys ${HOME}/.ssh/environment

# Binder keeps the ENTRYPOINT and passes its own Jupyter command as arguments; start.sh runs it
# in a login shell. For "docker exec" into a running container use: bash -lc '<command>'.
ENTRYPOINT ["/bin/bash", "/home/ubuntu/binder/start.sh"]
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser"]
