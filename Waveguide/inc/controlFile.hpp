#ifndef CONTROL_FILE_H
#define CONTROL_FILE_H
#include <string>
#include <jsoncpp/json/writer.h>
#include <vector>

/** Class that handles all the simulation parameters and saves them to a JSON file in the end */
class ControlFile
{
public:
  explicit ControlFile( const char* name );
  ControlFile(); // Use this if the control file should be loaded from another file
  ~ControlFile();

  /** Load an existing control file */
  void load( const std::string &fname );

  /** Save the parameters to file */
  void save() const;

  /** Add a field that is required */
  void addRequiredField( const char* reqField );

  /** Check that the file is valid. That is all the required fields are present */
  bool isValid() const; // Check that the required fields are present

  /** Get the top level JSON object */
  const Json::Value& get() const{ return *base; };

  /** Get the top level JSON object */
  Json::Value& get() { return *base; };

  /** Get a template for the filename. All files outputted should start with this */
  std::string getFnameTemplate() const { return fname; };

  /** Get the UID of the current run */
  unsigned int getUID() const { return uid; };
private:
  std::string fname;
  Json::Value *base;
  unsigned int uid;
  static const unsigned int uidDigits = 6;
  std::vector<std::string> requiredFields;
};

#endif
