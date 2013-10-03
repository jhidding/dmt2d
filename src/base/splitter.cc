#include "splitter.hh"

using namespace Misc;

Splitter::Splitter(std::string const &s, char delim_):
	data(s), delim(delim_) {}

Splitter::const_iterator Splitter::begin() const
{
	return const_iterator(*this, 0);
}

Splitter::const_iterator Splitter::end() const
{
	return const_iterator(*this);
}

Splitter::const_iterator::const_iterator(Splitter const &S_, size_t p):
	S(S_), p1(p) 
{
	if (p == std::string::npos)
		p2 = p1;
	else
		p2 = S.data.find(S.delim, p1);
}

Splitter::const_iterator &Splitter::const_iterator::operator++()
{
	if (p2 == S.data.length())
	{
		p1 = std::string::npos;
		return *this;
	}

	size_t a = S.data.find(S.delim, p2 + 1);
	if (a == std::string::npos) a = S.data.length();
	p1 = p2 + 1; p2 = a;

	return *this;
}

std::string Splitter::const_iterator::operator*() const
{
	if (p2 == std::string::npos)
		return S.data.substr(p1);

	return S.data.substr(p1, p2 - p1);
}

bool Splitter::const_iterator::operator==(const_iterator const &o) const
{
	return o.p1 == p1;
}

bool Splitter::const_iterator::operator!=(const_iterator const &o) const
{
	return o.p1 != p1;
}

