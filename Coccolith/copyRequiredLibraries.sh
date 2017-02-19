# Copies the required libraries to server
serverpath=/home/davidkl/privateLib/lib
serveraddress=davidkl@neutron.phys.ntnu.no

libPATH64="/usr/lib/x86_64-linux-gnu/"
libPATHLocal="/usr/local/lib/"
# Copy SFML
scp ${libPATH64}libsfml-* ${libPATH64}libfltk*  ${libPATHLocal}libpei.so ${libPATHLocal}libvisa.so ${serveraddress}:${serverpath}
