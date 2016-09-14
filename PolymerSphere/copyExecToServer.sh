source servers.sh
FILES="sphereScat.out sphere.msh"

scp ${FILES} davidkl@${IP}:${SERVER_PATH}
echo "Files copied to IP: ${IP}"
