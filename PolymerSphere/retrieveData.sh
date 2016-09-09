IP=129.241.184.56
FOLDER=/media/davidkl/stud/Master/OptiX/PolymerSphere/data
DESTINATION=data

scp davidkl@${IP}:${FOLDER}/*.bin ${DESTINATION}
scp davidkl@${IP}:${FOLDER}/*.json ${DESTINATION}
