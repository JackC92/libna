#ifndef NA_TYPE_TRAITS_NONESUCH_H
#define NA_TYPE_TRAITS_NONESUCH_H

namespace na
{
	struct nonesuch
	{
		nonesuch() = delete;
		~nonesuch() = delete;
		nonesuch(const nonesuch&) = delete;
		void operator=(const nonesuch&) = delete;
	};
}

#endif // !NA_TYPE_TRAITS_NONESUCH_H
