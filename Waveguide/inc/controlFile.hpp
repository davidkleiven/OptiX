#ifndef CONTROL_FILE_H
#define CONTROL_FILE_H
#include <string>
#include <jsoncpp/json/writer.h>
#include <vector>

class ControlFile
{
public:
  explicit ControlFile( const char* name );
  ControlFile(); // Use this if the control file should be loaded from another file
  ~ControlFile();
  void load( const std::string &fname );
  void save() const;
  void addRequiredField( const char* reqField );
  bool isValid() const; // Check that the required fields are present
  const Json::Value& get() const{ return *base; };
  Json::Value& get() { return *base; };
  std::string getFnameTemplate() const { return fname; };
  unsigned int getUID() const { return uid; };
private:
  std::string fname;
  Json::Value *base;
  unsigned int uid;
  static const unsigned int uidDigits = 6;
  std::vector<std::string> requiredFields;
};

#endif
