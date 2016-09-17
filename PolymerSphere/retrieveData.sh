source servers.sh
DESTINATION=data

key="$1"
if [ "$key" != "--uid" ]
then
  echo "Unknown key $key "
  exit 1
fi

uid="$2"
scp davidkl@${IP}:"${SERVER_PATH}/data/*$uid*" ${DESTINATION}
