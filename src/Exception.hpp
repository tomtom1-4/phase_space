#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <iostream>
#include <string>

class Exception : public std::exception {
  public:
    const char* what() const noexcept override {
      return "An error occured: ";
    }
};

#endif