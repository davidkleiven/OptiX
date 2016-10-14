#ifndef CONTROL_FILE_H
#define CONTROL_FILE_H
#include <string>
#include <jsoncpp/json/writer.h>

class ControlFile
{
public:
  explicit ControlFile( const char* name );
  ~ControlFile();
  void save() const;
  const Json::Value& get() const{ return *base; };
  Json::Value& get() { return *base; };
  std::string getFnameTemplate() const { return fname; };
private:
  std::string fname;
  Json::Value *base;
  unsigned int uid;
  static const unsigned int uidDigits = 6;
};

#endif
