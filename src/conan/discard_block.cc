#include "discard_block.h"

using namespace std;

void Conan::discard_block(istream &f)
{
	unsigned bs1, bs2;

	f.read(reinterpret_cast<char *>(&bs1), sizeof(unsigned));
	f.seekg(bs1, ios::cur);
	f.read(reinterpret_cast<char *>(&bs2), sizeof(unsigned));
	if (bs1 != bs2)
		throw "Could not discard a block";

/*	f.read(reinterpret_cast<char *>(&bs1), sizeof(unsigned));
	f.seekg(bs1, ios::cur);
	f.read(reinterpret_cast<char *>(&bs2), sizeof(unsigned));
	if (bs1 != bs2)
		throw "Could not discard a block";*/
}

