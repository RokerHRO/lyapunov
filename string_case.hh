#ifndef STRING_CASE_HH
#define STRING_CASE_HH

#include <string>
#include <string_view>

// replace all occurrences of search with replacement
inline
std::string replace_all(std::string orig, const std::string_view search, const std::string_view replacement)
{
	std::string::size_type pos = orig.find(search);
	while (pos < orig.size())
	{
		orig.replace(pos, search.length(), replacement);
		pos = orig.find(search, pos + replacement.length());
	}

	return orig;
}


constexpr
unsigned long long case_hash(const char* const s, const size_t len)
{
	auto uppercase = [](unsigned char c) -> unsigned { return c>='a' && c<='z' ? c-32 : c; };
	unsigned long long hash = 0;
	for (size_t i = 0; i < len; ++i)
	{
		hash = (hash * 67) + uppercase(s[i]) + 3;
	}
	return hash;
}

constexpr auto operator""_case(const char* const s, const size_t len)
{
	return case_hash(s, len);
}

constexpr auto case_hash(const std::string_view sv)
{
	return case_hash(sv.data(), sv.size());
}

#endif // STRING_CASE_HH
