source servers.sh
DESTINATION=data

scp davidkl@${IP}:${SERVER_PATH}/data/*.bin ${DESTINATION}
scp davidkl@${IP}:${SERVER_PATH}/data/*.json ${DESTINATION}
