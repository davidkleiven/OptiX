DDIR="dataPlane/MultInc30/WithEps"
PREFIX="ez"
PREFIX="${DDIR}/${PREFIX}"

# Convert all hdf5 files to png
echo "Converting hdf5 files to png..."
for file in ${PREFIX}*.h5;
do
  h5topng ${file} -Zc dkbluered
done;

# Convert png files to movie
echo "Converting pngs to animated gif..."
convert ${PREFIX}*.png ${PREFIX}.gif 
