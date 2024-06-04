#ifndef HelloCpp_H
#define HelloCpp_H

#include <string>
#include <vector>

class HelloCpp {
 private:
  std::string recipient;

 public:
  HelloCpp(std::string recipient = "World");
  void sayHello();
  void setRecipient(std::string recipient);
};

class VectorSource {
  private:
    std::vector<int> data;
  public:
    VectorSource(int size = 0);
    std::vector<int> getData();
    void set_size(int size);
};

#endif  // HelloCpp_H