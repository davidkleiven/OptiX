
void convertHDF5toPng( const string& h5file, const string& pngfile, const string& epsfile )
{
  string cmd("h5topng -S3 -Zc dkbluered -m -1.0 -M 1.0 -a yarg -A ");
  cmd += epsfile;
  cmd += " -o ";
  cmd += pngfile;
  cmd += " ";
  cmd += h5file;
  FILE* pipe = popen( cmd.c_str(), "r" );
  char buffer[128];
  string result("");
  if ( !pipe )
  {
    cerr << "Could not open pipe\n";
    return;
  }
  try
  {
    while( !feof(pipe) )
    {
      if( fgets( buffer, 128, pipe ) != NULL )
      {
        result += buffer;
      }
    }
  }
  catch(...)
  {
    pclose( pipe );
    return;
  }
  pclose(pipe);
}
