DDIR="dataPlane/MultInc30/WithEps"
PREFIX="ez"
PREFIX="${DDIR}/${PREFIX}"

# Convert all hdf5 files to png
echo "Converting hdf5 files to png..."
for file in ${PREFIX}*.h5;
do
  h5topng ${file} -Zc dkbluered
  FILENOEXT=${file%.*}
  PNGFILE="${FILENOEXT}.png"
  JPGFILE="${FILENOEXT}.jpg"
  convert ${PNGFILE} -resize 800x400 ${JPGFILE}
  rm ${PNGFILE}
done;

# Convert png files to movie
echo "Converting jpgs to animated gif..."
convert ${PREFIX}*.jpg ${PREFIX}.gif 
rm ${PREFIX}*.jpg
