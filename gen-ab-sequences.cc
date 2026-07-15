#include <iostream>
#include <string>
#include <deque>
#include <ranges>

std::deque<std::string> sequences;
unsigned base = 2;

std::string make_ab(unsigned long long u)
{
	std::string ab;
	while(u)
	{
		ab = char('A'+(u%base))  + ab;
		u = u/base;
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

std::string swap_ab(std::string ab)
{
	for (char& c : ab) { if (c=='A') c='#'; }
	for (char& c : ab) { if (c=='B') c='A'; }
	for (char& c : ab) { if (c=='#') c='B'; }
	return ab;
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


int main(int argc, char** argv)
{
	if (argc>1)
	{
		base = atoi(argv[1]);
		if (base<2 || base>6)
			base=2;
	}

	for(unsigned long long u=1;; ++u)
	{
		std::string ab = make_ab(u);
		if (ab.find('A') == std::string::npos)
			continue;
		if (ab.find('B') == std::string::npos)
			continue;
		if (is_repetition(ab))
			continue;

		const auto swapped = swap_ab(ab);
		if (std::ranges::find(sequences, swapped ) != sequences.end())
			continue;

		std::cout << ab << std::endl;
		sequences.emplace_back(std::move(ab));
	}
}
