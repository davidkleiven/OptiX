# Configuration script

OFNAME="searchPaths.sh"
HELP_MSG=$'Usage:./config.sh --baheaders <pathToHeaders> --help\n'
HELP_MSG+=$'--help: Print this message\n'
HELP_MSG+=$'--baheaders: Absolute path to the header file of the BornAgain library\n'
HELP_MSG+=$'--balib: Absolute path to the library files of the BornAgain library\n'

BA_HEADERS="./"
BA_LIB="./"
while [[ $# -gt 0 ]]
do
key="$1"
case $key in --baheaders)
  BA_HEADERS="$2"
  shift
  ;;
  --balib)
  BA_LIB="$2"
  shift
  ;;
  --help)
  echo "${HELP_MSG}"
  shift
  exit
  ;;
*)
  echo "Unknown option"
  ;;
esac
shift
done

# Append new paths to the search path file
echo "BA_HEADERS=${BA_HEADERS}" >> ${OFNAME}
echo "BA_LIB=${BA_LIB}" >> $OFNAME

echo "${OFNAME} is updated..."
