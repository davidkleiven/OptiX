source server.sh
#FILES="intensity2D_FD"

#for i in "${!FILES[@]}"
#do
#  FILES[$i]=${SERVER_PATH}/data/${FILES[$i]}
#done
scp davidkl@${IP}:${SERVER_PATH}/data/*$1* data/
