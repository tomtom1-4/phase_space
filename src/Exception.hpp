#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <iostream>
#include <string>

class Exception : public std::exception {
  public:
    std::string message;

    const char* what() const noexcept override {
      std::cout << "Exception: " << std::endl;
      return message.c_str();
    }

    Exception(const std::string& message) : message(message) { }
};

#endif