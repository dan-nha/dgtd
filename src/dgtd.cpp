#include <iostream>

void stream_welcome_message() {
  std::cout << "\033[31m"
    << "---------------------------------------\n"
    << "Welcome to DGTD \n\n"
    << "This is a one-dimensional \n"
    << "discontinuous Galerkin time-domain \n"
    << "(DGTD) solver.\n"
    << "---------------------------------------\n\033[0m"
    << std::endl;
}

int main(int argc, char* argv[]) {
  stream_welcome_message();
}

