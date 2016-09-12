source servers.sh
DESTINATION=data

scp davidkl@${IP}:"${SERVER_PATH}/data/*.json ${SERVER_PATH}/data/*.bin ${SERVER_PATH}/data/*.pp ${SERVER_PATH}/data/*.h5" ${DESTINATION}
