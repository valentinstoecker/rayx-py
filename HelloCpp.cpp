#include "HelloCpp.h"

#include <iostream>

HelloCpp::HelloCpp(std::string recipient) : recipient(recipient) {}

void HelloCpp::sayHello() {
  std::cout << "Hello, " << recipient << "!" << std::endl;
}

void HelloCpp::setRecipient(std::string recipient) {
  this->recipient = recipient;
}

VectorSource::VectorSource(int size) {
  set_size(size);
}

std::vector<int> VectorSource::getData() {
  return data;
}

void VectorSource::set_size(int size) {
  data.resize(size);
  for (int i = 0; i < size; i++) {
    data[i] = i;
  }
}
