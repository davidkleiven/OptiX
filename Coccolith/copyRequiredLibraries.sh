# Copies the required libraries to server
serverpath=/home/davidkl/privateLib/lib
serveraddress=davidkl@neutron.phys.ntnu.no

#libPATH64="/usr/lib/x86_64-linux-gnu/"
libPATH64="/usr/lib64/"
libPATHLocal="/usr/local/lib/"
H5Path="/usr/lib64/"

# Copy SFML
scp ${libPATH64}libsfml-* ${libPATH64}libfltk*  ${libPATHLocal}libpei.so ${libPATHLocal}libvisa.so \
  ${libPATH64}/libopenblas.so* ${H5Path}libhdf5* ${serveraddress}:${serverpath}
