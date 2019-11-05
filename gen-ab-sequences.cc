#include <iostream>
#include <string>
#include <deque>

std::deque<std::string> sequences;

std::string make_ab(unsigned long long u)
{
	std::string ab;
	while(u)
	{
		ab = (u&1 ? 'A' : 'B') + ab;
		u = u/2;
	}
	return ab;
}


std::string operator*(const std::string& o, unsigned count)
{
	std::string ret;
	for(unsigned u=0; u<count; ++u)
	{
		ret += o;
	}
	return ret;
}


// return true if 'ab' is a repetition of an already found sequence
bool is_repetition(const std::string& ab)
{
	const unsigned max_len = ab.size() / 2;
	for(const auto& s : sequences)
	{
		if(s.size() > max_len)
			return false;
		
		if(ab.size() % s.size())  // ab cannot be a repetition of s
			continue;
		
		const unsigned rep = ab.size() / s.size();
		
		if( s * rep == ab)
		{
			std::cerr << "\"" << ab << "\" = \"" << s << "\" * " << rep << ". SKIP!\n";
			return true;
		}
	}
	return false;
}


int main()
{
	for(unsigned long long u=1;; ++u)
	{
		std::string ab = make_ab(u);
		if(ab.find('B') != std::string::npos && !is_repetition(ab))
		{
			std::cout << ab << std::endl;
			sequences.emplace_back(std::move(ab));
		}
	}
}
