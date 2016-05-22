#ifndef DATA_TO_FILE_H
#define DATA_TO_FILE_H
#include <string>
#include <vector>
#include "meep.hpp"
typedef std::vector<double> dvec;
void monitorToFile(const std::string &fname, const dvec &time, const dvec &fieldVal, const meep::vec &pos);
#endif
