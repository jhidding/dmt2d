#include <typeinfo>
#include <utility>
#include <complex>

#include "mvector.h"
#include "mtypeid.h"

using namespace std;
using namespace Conan;

TypeIdMap	*TypeIdMap::singleton = NULL;
string	const	TypeIdMap::default_name = "unknown";

template<typename T>
inline pair<string, string> makeId(string const &id)
{
	return pair<string,string>(typeid(T).name(), id);
}

pair<string, string> prefab[] = 
{
	makeId<int>("int32"),
	makeId<unsigned>("uint32"),
	makeId<float>("float32"),
	makeId<double>("float64"),
	makeId<complex<float> >("complex64"),
	makeId<complex<double> >("complex128"),
	makeId<mVector<int, 3>>("3i4"),
	makeId<mVector<float, 3>>("3f4"),
	makeId<mVector<float, 2>>("2f4"),
	makeId<mVector<double, 3>>("3f8"),
	makeId<mVector<double, 2>>("2f8"),
};

TypeIdMap::TypeIdMap()
	: map<string,string>(prefab, prefab + 11)
{
}

string const &TypeIdMap::operator[](string const &idx)
{
	if (count(idx) == 1)
		return map<string,string>::operator[](idx);
	else
		return default_name;
}

TypeIdMap &TypeIdMap::instance()
{
	if (singleton == NULL)
		singleton = new TypeIdMap;

	return *singleton;
}

/*mAddTypeId<int>("int");
mAddTypeId<unsigned>("unsigned");
mAddTypeId<float>("float");
mAddTypeId<double>("double");
mAddTypeId<complex<double> >("complex<double>");
*/

