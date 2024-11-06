#include "headers/wave.h"

#include <limits>

namespace
{

	/**
	 * Return distribution * log(distribution).
	 */
	std::vector<double>
	get_plogp(const std::vector<double> &distribution) noexcept
	{
		std::vector<double> plogp;
		for (unsigned i = 0; i < distribution.size(); i++)
		{
			plogp.push_back(distribution[i] * log(distribution[i]));
		}
		return plogp;
	}

} // namespace

Wave::Wave(unsigned height, unsigned width,
		   const std::vector<double> &patterns_frequencies) noexcept
	: patterns_frequencies(patterns_frequencies),
	  plogp_patterns_frequencies(get_plogp(patterns_frequencies)),
	  is_impossible(false), nb_patterns(patterns_frequencies.size()),
	  data(width * height, nb_patterns, 1), width(width), height(height),
	  size(height * width)
{
	// Initialize the memoisation of entropy.
	double base_entropy = 0;
	double base_s = 0;
	for (unsigned i = 0; i < nb_patterns; i++)
	{
		base_entropy += plogp_patterns_frequencies[i];
		base_s += patterns_frequencies[i];
	}
	double log_base_s = log(base_s);
	double entropy_base = log_base_s - base_entropy / base_s;
	memoisation.plogp_sum = std::vector<double>(width * height, base_entropy);
	memoisation.sum = std::vector<double>(width * height, base_s);
	memoisation.log_sum = std::vector<double>(width * height, log_base_s);
	memoisation.nb_patterns = std::vector<unsigned>(width * height, static_cast<unsigned>(nb_patterns));
	memoisation.entropy = std::vector<double>(width * height, entropy_base);
	lowest_entropy = entropy_base;
	is_index_on_lowest = std::vector<bool>(width * height, false);
	lowest_indexes.reserve(width * height);
}

void Wave::set(unsigned index, unsigned pattern, bool value) noexcept
{
	bool old_value = data.get(index, pattern);
	// If the value isn't changed, nothing needs to be done.
	if (old_value == value)
	{
		return;
	}
	// Otherwise, the memoisation should be updated.
	data.get(index, pattern) = value;
	memoisation.plogp_sum[index] -= plogp_patterns_frequencies[pattern];
	memoisation.sum[index] -= patterns_frequencies[pattern];
	memoisation.log_sum[index] = log(memoisation.sum[index]);
	memoisation.nb_patterns[index]--;
	memoisation.entropy[index] =
		memoisation.log_sum[index] -
		memoisation.plogp_sum[index] / memoisation.sum[index];
	// If there is no patterns possible in the cell, then there is a
	// contradiction.
	if (memoisation.nb_patterns[index] == 0)
	{
		is_impossible = true;
	}
	// If the new entropy is equal to lowest, we add the index to lowest_indexes
	if (!is_index_on_lowest[index] && memoisation.entropy[index] == lowest_entropy)
	{
		lowest_indexes.push_back(index);
		is_index_on_lowest[index] = true;
	}
	// If the new entropy is lower than lowest, now this entropy value is the new
	// lowest and we reset is_index_on_lowest and lowest_indexes
	else if (memoisation.entropy[index] < lowest_entropy)
	{
		lowest_entropy = memoisation.entropy[index];
		lowest_indexes.clear();
		lowest_indexes.push_back(index);
		is_index_on_lowest = std::vector<bool>(size, false);
		is_index_on_lowest[index] = true;
	}
}

int Wave::get_min_entropy(std::minstd_rand &gen) noexcept
{
	if (is_impossible)
	{
		return -2;
	}

	int argmin = -1;

	// We check if there are lowest indexes that we know.
	// If not, then we'll find the new ones (the original way)
	if (lowest_indexes.empty())
	{
		lowest_entropy = std::numeric_limits<double>::infinity();
		for (unsigned i = 0; i < size; i++)
		{
			double nb_patterns_local = memoisation.nb_patterns[i];
			if (nb_patterns_local == 1)
			{
				continue;
			}
			// We take the memoised entropy.
			double entropy = memoisation.entropy[i];

			// If the new entropy is equal to lowest, we add the index to lowest_indexes
			if (entropy == lowest_entropy)
			{
				lowest_indexes.push_back(i);
				is_index_on_lowest[i] = true;
			}
			// If the new entropy is lower than lowest, now this entropy value is the new
			// lowest and we reset is_index_on_lowest and lowest_indexes
			else if (entropy < lowest_entropy)
			{
				lowest_entropy = entropy;
				lowest_indexes.clear();
				lowest_indexes.push_back(i);
				is_index_on_lowest = std::vector<bool>(size, false);
				is_index_on_lowest[i] = true;
			}
		}
	}

	// We pick one of the lowest indexes randomly if lowest indexes is not empty
	if (!lowest_indexes.empty())
	{
		unsigned rand = gen();
		unsigned selecte_index = rand % lowest_indexes.size();
		argmin = lowest_indexes[selecte_index];
		lowest_indexes.erase(lowest_indexes.begin() + selecte_index);
	}

	return argmin;
}
