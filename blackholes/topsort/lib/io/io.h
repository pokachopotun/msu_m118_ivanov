#pragma once

#include <common.h>

#include <iostream>
#include <fstream>
#include <string>

TGraph ReadGraphFromBinaryStream(std::istream& file);
TGraph ReadGraphFromBinaryFile(const std::string& filename);

TGraph ReadGraphFromTextStream(std::istream& stream);
TGraph ReadGraphFromTextFile(const std::string& filename);
