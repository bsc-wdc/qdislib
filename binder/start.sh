#!/bin/bash
# COMPSs starts its worker with "ssh -p 22 localhost". Binder runs as a non-root user, so start a
# private sshd on 127.0.0.1:2222 (the ssh wrapper in the Dockerfile redirects COMPSs there).
/usr/sbin/sshd -f "$HOME/.sshd/sshd_config" -E "$HOME/.sshd/sshd.log" || echo "warning: sshd did not start" >&2
# Run the command in a login shell so /etc/profile.d/compss.sh sets up the COMPSs environment.
exec /bin/bash -lc 'exec "$@"' _ "$@"
