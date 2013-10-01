#pragma once
#include <string>
#include <map>
#include <typeinfo>

namespace Conan
{

class TypeIdMap: public std::map<std::string, std::string>
{
	static TypeIdMap	*singleton;
	static std::string const default_name;

	public:
		static TypeIdMap &instance();

		std::string const &operator[](std::string const &idx);

	private:
		TypeIdMap();
};

template <typename T>
std::string mTypeId()
{
	return TypeIdMap::instance()[typeid(T).name()];
}

template <typename T>
void mAddTypeId(std::string const &id)
{
	TypeIdMap::instance().insert(std::pair<std::string,std::string>(
				typeid(T).name(),
				id));
}

}

