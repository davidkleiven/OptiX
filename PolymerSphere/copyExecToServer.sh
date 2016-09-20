source servers.sh
FILES="sphere.msh cube.msh"

executable="$1"
FILES="${FILES} $executable"
if [ -z $executable ]
then
  echo "No executable sepcified..."
  exit 1
fi

scp ${FILES} davidkl@${IP}:${SERVER_PATH}
echo "Files copied to IP: ${IP}"
