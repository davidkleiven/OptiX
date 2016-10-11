source server.sh

EXECUTABLE=$1
scp ${EXECUTABLE} davidkl@${IP}:${SERVER_PATH}/
