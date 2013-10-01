#pragma once
#include <string>

extern std::string default_filename(std::string const &q);
extern std::string adhesion_filename(float b);
extern std::string timed_filename(std::string const &q, float b = -1.0);

