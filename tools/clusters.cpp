// new cluster c file
//    /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/tools

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

int main()
{ 
  std::ifstream ifile("/uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/tools")

    while (std::getline(ifile, line)) // read the current line
      {
	std::istringstream iss{line}; // construct a string stream from line

        // read the tokens from current line separated by comma
	std::vector<std::string> tokens; // here we store the tokens
	std::string token; // current token
        while (std::getline(iss, token, ','))
	  {
            tokens.push_back(token); // add the token to the vector
	  }
	std::cout << "Tokenized line: ";
        for (const auto& elem : tokens)
	  std::cout << "[" << elem << "]";
	std::cout << std::endl;
      }
}
