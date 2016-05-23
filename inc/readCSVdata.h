#ifndef READ_CSV_DATA_H
#define READ_CSV_DATA_H
#include <string>
#include <vector>
class ReadCSVData
{
public:
  ReadCSVData():nCol(1){};
  void read(const std::string &fname, unsigned int nColumns);
  double get(unsigned int row, unsigned int col) const;
  unsigned int numPoints() const { return data.size()/nCol; };
private:
  std::vector<double> data;
  unsigned int nCol;
};
#endif
