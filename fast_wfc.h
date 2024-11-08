#ifndef FAST_WFC_H
#define FAST_WFC_H

#include "assert.h"

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <array>
#include <random>

constexpr int directions_x[4] = {0, -1, 1, 0};
constexpr int directions_y[4] = {-1, 0, 0, 1};

constexpr unsigned get_opposite_direction(unsigned direction) noexcept
{
	return 3 - direction;
}

template <typename T>
class Array2D
{

public:
	std::size_t height;
	std::size_t width;

	std::vector<T> data;

	Array2D(std::size_t _height, std::size_t _width) noexcept
		: height(_height), width(_width), data(_width * _height) {}

	Array2D(std::size_t _height, std::size_t _width, T value) noexcept
		: height(_height), width(_width), data(_width * _height, value) {}

	const T &get(std::size_t i, std::size_t j) const noexcept
	{
		assert(i < height && j < width);
		return data[j + i * width];
	}

	T &get(std::size_t i, std::size_t j) noexcept
	{
		assert(i < height && j < width);
		return data[j + i * width];
	}

	Array2D<T> reflected() const noexcept
	{
		Array2D<T> result = Array2D<T>(width, height);
		for (std::size_t y = 0; y < height; y++)
		{
			for (std::size_t x = 0; x < width; x++)
			{
				result.get(y, x) = get(y, width - 1 - x);
			}
		}
		return result;
	}

	Array2D<T> rotated() const noexcept
	{
		Array2D<T> result = Array2D<T>(width, height);
		for (std::size_t y = 0; y < width; y++)
		{
			for (std::size_t x = 0; x < height; x++)
			{
				result.get(y, x) = get(x, width - 1 - y);
			}
		}
		return result;
	}

	Array2D<T> get_sub_array(std::size_t y, std::size_t x, std::size_t sub_width,
							 std::size_t sub_height) const noexcept
	{
		Array2D<T> sub_array_2d = Array2D<T>(sub_width, sub_height);
		for (std::size_t ki = 0; ki < sub_height; ki++)
		{
			for (std::size_t kj = 0; kj < sub_width; kj++)
			{
				sub_array_2d.get(ki, kj) = get((y + ki) % height, (x + kj) % width);
			}
		}
		return sub_array_2d;
	}

	bool operator==(const Array2D<T> &a) const noexcept
	{
		if (height != a.height || width != a.width)
		{
			return false;
		}

		for (std::size_t i = 0; i < data.size(); i++)
		{
			if (a.data[i] != data[i])
			{
				return false;
			}
		}
		return true;
	}
};

namespace std
{
	template <typename T>
	class hash<Array2D<T>>
	{
	public:
		std::size_t operator()(const Array2D<T> &a) const noexcept
		{
			std::size_t seed = a.data.size();
			for (const T &i : a.data)
			{
				seed ^= hash<T>()(i) + (std::size_t)0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};
}

template <typename T>
class Array3D
{

public:
	std::size_t height;
	std::size_t width;
	std::size_t depth;

	std::vector<T> data;

	Array3D(std::size_t _height, std::size_t _width, std::size_t _depth) noexcept
		: height(_height), width(_width), depth(_depth),
		  data(_width * _height * _depth) {}

	Array3D(std::size_t _height, std::size_t _width, std::size_t _depth,
			T value) noexcept
		: height(_height), width(_width), depth(_depth),
		  data(_width * _height * _depth, value) {}

	const T &get(std::size_t i, std::size_t j, std::size_t k) const noexcept
	{
		assert(i < height && j < width && k < depth);
		return data[i * width * depth + j * depth + k];
	}

	T &get(std::size_t i, std::size_t j, std::size_t k) noexcept
	{
		return data[i * width * depth + j * depth + k];
	}

	bool operator==(const Array3D &a) const noexcept
	{
		if (height != a.height || width != a.width || depth != a.depth)
		{
			return false;
		}

		for (std::size_t i = 0; i < data.size(); i++)
		{
			if (a.data[i] != data[i])
			{
				return false;
			}
		}
		return true;
	}
};

struct OverlappingWFCOptions
{
	bool periodic_input;
	bool periodic_output;
	unsigned out_height;
	unsigned out_width;
	unsigned symmetry;
	bool ground;
	unsigned pattern_size;

	unsigned get_wave_height() const noexcept
	{
		return periodic_output ? out_height : out_height - pattern_size + 1;
	}

	unsigned get_wave_width() const noexcept
	{
		return periodic_output ? out_width : out_width - pattern_size + 1;
	}
};

struct EntropyMemoisation
{
	std::vector<double> plogp_sum;
	std::vector<double> sum;
	std::vector<double> log_sum;
	std::vector<unsigned> nb_patterns;
	std::vector<double> entropy;
};

class Wave
{
private:
	const std::vector<double> patterns_frequencies;

	const std::vector<double> plogp_patterns_frequencies;

	const double min_abs_half_plogp;

	EntropyMemoisation memoisation;

	bool is_impossible;

	const size_t nb_patterns;

	Array2D<uint8_t> data;

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

	double get_min_abs_half(const std::vector<double> &v) noexcept
	{
		double min_abs_half = std::numeric_limits<double>::infinity();
		for (unsigned i = 0; i < v.size(); i++)
		{
			min_abs_half = std::min(min_abs_half, std::abs(v[i] / 2.0));
		}
		return min_abs_half;
	}

public:
	const unsigned width;
	const unsigned height;
	const unsigned size;

	Wave(
		unsigned height,
		unsigned width,
		const std::vector<double> &patterns_frequencies) noexcept
		: patterns_frequencies(patterns_frequencies),
		  plogp_patterns_frequencies(get_plogp(patterns_frequencies)),
		  min_abs_half_plogp(get_min_abs_half(plogp_patterns_frequencies)),
		  is_impossible(false),
		  nb_patterns(patterns_frequencies.size()),
		  data(width * height, nb_patterns, 1),
		  width(width),
		  height(height),
		  size(height * width)
	{
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
	}

	bool get(unsigned index, unsigned pattern) const noexcept
	{
		return data.get(index, pattern);
	}

	bool get(unsigned i, unsigned j, unsigned pattern) const noexcept
	{
		return get(i * width + j, pattern);
	}

	void set(unsigned index, unsigned pattern, bool value) noexcept
	{
		bool old_value = data.get(index, pattern);
		if (old_value == value)
		{
			return;
		}
		data.get(index, pattern) = value;
		memoisation.plogp_sum[index] -= plogp_patterns_frequencies[pattern];
		memoisation.sum[index] -= patterns_frequencies[pattern];
		memoisation.log_sum[index] = log(memoisation.sum[index]);
		memoisation.nb_patterns[index]--;
		memoisation.entropy[index] =
			memoisation.log_sum[index] -
			memoisation.plogp_sum[index] / memoisation.sum[index];
		if (memoisation.nb_patterns[index] == 0)
		{
			is_impossible = true;
		}
	}

	void set(unsigned i, unsigned j, unsigned pattern, bool value) noexcept
	{
		set(i * width + j, pattern, value);
	}

	int get_min_entropy(std::minstd_rand &gen) const noexcept
	{
		if (is_impossible)
		{
			return -2;
		}

		std::uniform_real_distribution<> dis(0, min_abs_half_plogp);

		double min = std::numeric_limits<double>::infinity();
		int argmin = -1;

		for (unsigned i = 0; i < size; i++)
		{
			double nb_patterns_local = memoisation.nb_patterns[i];
			if (nb_patterns_local == 1)
			{
				continue;
			}

			double entropy = memoisation.entropy[i];

			if (entropy <= min)
			{
				double noise = dis(gen);
				if (entropy + noise < min)
				{
					min = entropy + noise;
					argmin = i;
				}
			}
		}

		return argmin;
	}
};

class Propagator
{
public:
	using PropagatorState = std::vector<std::array<std::vector<unsigned>, 4>>;

private:
	const std::size_t patterns_size;

	PropagatorState propagator_state;

	const unsigned wave_width;
	const unsigned wave_height;

	const bool periodic_output;

	std::vector<std::tuple<unsigned, unsigned, unsigned>> propagating;

	Array3D<std::array<int, 4>> compatible;

	void init_compatible() noexcept
	{
		std::array<int, 4> value;
		for (unsigned y = 0; y < wave_height; y++)
		{
			for (unsigned x = 0; x < wave_width; x++)
			{
				for (unsigned pattern = 0; pattern < patterns_size; pattern++)
				{
					for (int direction = 0; direction < 4; direction++)
					{
						value[direction] = static_cast<unsigned>(
							propagator_state[pattern][get_opposite_direction(direction)].size());
					}
					compatible.get(y, x, pattern) = value;
				}
			}
		}
	}

public:
	Propagator(unsigned wave_height, unsigned wave_width, bool periodic_output,
			   PropagatorState propagator_state) noexcept
		: patterns_size(propagator_state.size()),
		  propagator_state(propagator_state), wave_width(wave_width),
		  wave_height(wave_height), periodic_output(periodic_output),
		  compatible(wave_height, wave_width, patterns_size)
	{
		init_compatible();
	}

	void add_to_propagator(unsigned y, unsigned x, unsigned pattern) noexcept
	{
		std::array<int, 4> temp = {};
		compatible.get(y, x, pattern) = temp;
		propagating.emplace_back(y, x, pattern);
	}

	void propagate(Wave &wave) noexcept
	{
		while (propagating.size() != 0)
		{
			unsigned y1, x1, pattern;
			std::tie(y1, x1, pattern) = propagating.back();
			propagating.pop_back();

			for (unsigned direction = 0; direction < 4; direction++)
			{
				int dx = directions_x[direction];
				int dy = directions_y[direction];
				int x2, y2;
				if (periodic_output)
				{
					x2 = ((int)x1 + dx + (int)wave.width) % wave.width;
					y2 = ((int)y1 + dy + (int)wave.height) % wave.height;
				}
				else
				{
					x2 = x1 + dx;
					y2 = y1 + dy;
					if (x2 < 0 || x2 >= (int)wave.width)
					{
						continue;
					}
					if (y2 < 0 || y2 >= (int)wave.height)
					{
						continue;
					}
				}

				unsigned i2 = x2 + y2 * wave.width;
				const std::vector<unsigned> &patterns =
					propagator_state[pattern][direction];

				for (auto it = patterns.begin(), it_end = patterns.end(); it < it_end;
					 ++it)
				{
					std::array<int, 4> &value = compatible.get(y2, x2, *it);
					value[direction]--;

					if (value[direction] == 0)
					{
						add_to_propagator(y2, x2, *it);
						wave.set(i2, *it, false);
					}
				}
			}
		}
	}
};

class WFC
{
private:
	std::minstd_rand gen;

	const std::vector<double> patterns_frequencies;

	Wave wave;

	const size_t nb_patterns;

	Propagator propagator;

	Array2D<unsigned> *wave_to_output() const noexcept
	{
		Array2D<unsigned> *output_patterns = new Array2D<unsigned>(wave.height, wave.width);
		for (unsigned i = 0; i < wave.size; i++)
		{
			for (unsigned k = 0; k < nb_patterns; k++)
			{
				if (wave.get(i, k))
				{
					output_patterns->data[i] = k;
				}
			}
		}
		return output_patterns;
	}

	std::vector<double> &normalize(std::vector<double> &v)
	{
		double sum_weights = 0.0;
		for (double weight : v)
		{
			sum_weights += weight;
		}

		double inv_sum_weights = 1.0 / sum_weights;
		for (double &weight : v)
		{
			weight *= inv_sum_weights;
		}

		return v;
	}

public:
	WFC(
		bool periodic_output,
		unsigned seed,
		std::vector<double> patterns_frequencies,
		Propagator::PropagatorState propagator,
		unsigned wave_height,
		unsigned wave_width) noexcept
		: gen(seed),
		  patterns_frequencies(normalize(patterns_frequencies)),
		  wave(wave_height, wave_width, patterns_frequencies),
		  nb_patterns(propagator.size()),
		  propagator(wave.height, wave.width, periodic_output, propagator)
	{
	}

	Array2D<unsigned> *run() noexcept
	{
		while (true)
		{
			ObserveStatus result = observe();

			if (result == FAILURE)
			{
				return nullptr;
			}
			else if (result == SUCCESS)
			{
				return wave_to_output();
			}

			propagator.propagate(wave);
		}
	}

	enum ObserveStatus
	{
		SUCCESS,
		FAILURE,
		TO_CONTINUE
	};

	ObserveStatus observe() noexcept
	{
		int argmin = wave.get_min_entropy(gen);

		if (argmin == -2)
		{
			return FAILURE;
		}

		if (argmin == -1)
		{
			wave_to_output();
			return SUCCESS;
		}

		double s = 0;
		for (unsigned k = 0; k < nb_patterns; k++)
		{
			s += wave.get(argmin, k) ? patterns_frequencies[k] : 0;
		}

		std::uniform_real_distribution<> dis(0, s);
		double random_value = dis(gen);
		size_t chosen_value = nb_patterns - 1;

		for (unsigned k = 0; k < nb_patterns; k++)
		{
			random_value -= wave.get(argmin, k) ? patterns_frequencies[k] : 0;
			if (random_value <= 0)
			{
				chosen_value = k;
				break;
			}
		}

		for (unsigned k = 0; k < nb_patterns; k++)
		{
			if (wave.get(argmin, k) != (k == chosen_value))
			{
				propagator.add_to_propagator(argmin / wave.width, argmin % wave.width,
											 k);
				wave.set(argmin, k, false);
			}
		}

		return TO_CONTINUE;
	}

	void propagate() noexcept { propagator.propagate(wave); }

	void remove_wave_pattern(unsigned i, unsigned j, unsigned pattern) noexcept
	{
		if (wave.get(i, j, pattern))
		{
			wave.set(i, j, pattern, false);
			propagator.add_to_propagator(i, j, pattern);
		}
	}
};

class OverlappingWFC
{

private:
	Array2D<uint32_t> input;

	OverlappingWFCOptions options;

	std::vector<Array2D<uint32_t>> patterns;

	WFC wfc;

	OverlappingWFC(
		const Array2D<uint32_t> &_input,
		const OverlappingWFCOptions &_options,
		const unsigned &_seed,
		const std::pair<std::vector<Array2D<uint32_t>>, std::vector<double>> &_patterns,
		const Propagator::PropagatorState &_propagator) noexcept
		: input(_input),
		  options(_options),
		  patterns(_patterns.first),
		  wfc(_options.periodic_output, _seed, _patterns.second, _propagator, _options.get_wave_height(), _options.get_wave_width())
	{
		if (_options.ground)
		{
			init_ground(wfc, _input, _patterns.first, _options);
		}
	}

	OverlappingWFC(
		const Array2D<uint32_t> &_input,
		const OverlappingWFCOptions &_options,
		const unsigned &_seed,
		const std::pair<std::vector<Array2D<uint32_t>>, std::vector<double>> &_patterns) noexcept
		: OverlappingWFC(
			  _input, _options, _seed, _patterns,
			  generate_compatible(_patterns.first)) {}

	void init_ground(
		WFC &_wfc,
		const Array2D<uint32_t> &_input,
		const std::vector<Array2D<uint32_t>> &_patterns,
		const OverlappingWFCOptions &_options) noexcept
	{
		unsigned ground_pattern_id = get_ground_pattern_id(_input, _patterns, _options);

		for (unsigned j = 0; j < _options.get_wave_width(); j++)
		{
			set_pattern(ground_pattern_id, _options.get_wave_height() - 1, j);
		}

		for (unsigned i = 0; i < _options.get_wave_height() - 1; i++)
		{
			for (unsigned j = 0; j < _options.get_wave_width(); j++)
			{
				wfc.remove_wave_pattern(i, j, ground_pattern_id);
			}
		}

		_wfc.propagate();
	}

	static unsigned
	get_ground_pattern_id(
		const Array2D<uint32_t> &_input,
		const std::vector<Array2D<uint32_t>> &_patterns,
		const OverlappingWFCOptions &_options) noexcept
	{
		std::size_t y = _input.height - 1;
		std::size_t x = _input.width / 2;
		unsigned size = _options.pattern_size;
		Array2D<uint32_t> ground_pattern = _input.get_sub_array(y, x, size, size);

		for (unsigned i = 0; i < _patterns.size(); i++)
		{
			if (ground_pattern == _patterns[i])
			{
				return i;
			}
		}

		assert(false);
		return 0;
	}

	static std::pair<std::vector<Array2D<uint32_t>>, std::vector<double>>
	get_patterns(
		const Array2D<uint32_t> &_input,
		const OverlappingWFCOptions &_options) noexcept
	{
		std::vector<Array2D<uint32_t>> _patterns;
		std::unordered_map<Array2D<uint32_t>, unsigned> patterns_id;

		std::vector<double> patterns_weight;

		unsigned size = _options.pattern_size;
		std::vector<Array2D<uint32_t>> symmetries(8, Array2D<uint32_t>(size, size));

		unsigned max_i = _options.periodic_input ? _input.height : _input.height - size + 1;
		unsigned max_j = _options.periodic_input ? _input.width : _input.width - size + 1;

		for (unsigned i = 0; i < max_i; i++)
		{
			for (unsigned j = 0; j < max_j; j++)
			{
				symmetries[0].data = _input.get_sub_array(i, j, size, size).data;
				symmetries[1].data = symmetries[0].reflected().data;
				symmetries[2].data = symmetries[0].rotated().data;
				symmetries[3].data = symmetries[2].reflected().data;
				symmetries[4].data = symmetries[2].rotated().data;
				symmetries[5].data = symmetries[4].reflected().data;
				symmetries[6].data = symmetries[4].rotated().data;
				symmetries[7].data = symmetries[6].reflected().data;

				for (unsigned k = 0; k < _options.symmetry; k++)
				{
					auto res = patterns_id.insert(
						std::make_pair(symmetries[k], _patterns.size()));

					if (!res.second)
					{
						patterns_weight[res.first->second] += 1;
					}
					else
					{
						_patterns.push_back(symmetries[k]);
						patterns_weight.push_back(1);
					}
				}
			}
		}
		return {_patterns, patterns_weight};
	}

	static bool agrees(const Array2D<uint32_t> &pattern1, const Array2D<uint32_t> &pattern2, int dy, int dx) noexcept
	{
		unsigned xmin = dx < 0 ? 0 : dx;
		unsigned xmax = dx < 0 ? dx + pattern2.width : pattern1.width;
		unsigned ymin = dy < 0 ? 0 : dy;
		unsigned ymax = dy < 0 ? dy + pattern2.height : pattern1.width;

		for (unsigned y = ymin; y < ymax; y++)
		{
			for (unsigned x = xmin; x < xmax; x++)
			{
				if (pattern1.get(y, x) != pattern2.get(y - dy, x - dx))
				{
					return false;
				}
			}
		}
		return true;
	}

	static std::vector<std::array<std::vector<unsigned>, 4>>
	generate_compatible(const std::vector<Array2D<uint32_t>> &_patterns) noexcept
	{
		std::vector<std::array<std::vector<unsigned>, 4>> compatible =
			std::vector<std::array<std::vector<unsigned>, 4>>(_patterns.size());

		for (unsigned p1 = 0; p1 < _patterns.size(); p1++)
		{
			for (unsigned d = 0; d < 4; d++)
			{
				for (unsigned p2 = 0; p2 < _patterns.size(); p2++)
				{
					if (agrees(_patterns[p1], _patterns[p2], directions_y[d], directions_x[d]))
					{
						compatible[p1][d].push_back(p2);
					}
				}
			}
		}

		return compatible;
	}

	Array2D<uint32_t> *to_image(const Array2D<unsigned> *output_patterns) const noexcept
	{
		Array2D<uint32_t> &output = *new Array2D<uint32_t>(options.out_height, options.out_width);

		unsigned height = options.get_wave_height();
		unsigned width = options.get_wave_width();
		if (options.periodic_output)
		{
			for (unsigned y = 0; y < height; y++)
			{
				for (unsigned x = 0; x < width; x++)
				{
					output.get(y, x) = patterns[output_patterns->get(y, x)].get(0, 0);
				}
			}
		}
		else
		{
			for (unsigned y = 0; y < height; y++)
			{
				for (unsigned x = 0; x < width; x++)
				{
					output.get(y, x) = patterns[output_patterns->get(y, x)].get(0, 0);
				}
			}
			for (unsigned y = 0; y < height; y++)
			{
				const Array2D<unsigned> &pattern = patterns[output_patterns->get(y, width - 1)];
				for (unsigned dx = 1; dx < options.pattern_size; dx++)
				{
					output.get(y, width - 1 + dx) = pattern.get(0, dx);
				}
			}
			for (unsigned x = 0; x < width; x++)
			{
				const Array2D<unsigned> &pattern = patterns[output_patterns->get(height - 1, x)];
				for (unsigned dy = 1; dy < options.pattern_size; dy++)
				{
					output.get(height - 1 + dy, x) = pattern.get(dy, 0);
				}
			}
			const Array2D<unsigned> &pattern = patterns[output_patterns->get(height - 1, width - 1)];
			for (unsigned dy = 1; dy < options.pattern_size; dy++)
			{
				for (unsigned dx = 1; dx < options.pattern_size; dx++)
				{
					output.get(height - 1 + dy, width - 1 + dx) = pattern.get(dy, dx);
				}
			}
		}

		return &output;
	}

	int get_pattern_id(const Array2D<unsigned> &_pattern)
	{
		std::vector<Array2D<unsigned>>::iterator it = std::find(patterns.begin(), patterns.end(), _pattern);
		if (it != patterns.end())
			return patterns.end() - it;
		else
			return -1;
	}

	void set_pattern(unsigned pattern_id, unsigned i, unsigned j) noexcept
	{
		for (unsigned p = 0; p < patterns.size(); p++)
		{
			if (pattern_id != p)
			{
				wfc.remove_wave_pattern(i, j, p);
			}
		}
	}

public:
	OverlappingWFC(
		const Array2D<uint32_t> &_input,
		const OverlappingWFCOptions &_options,
		int _seed) noexcept
		: OverlappingWFC(_input, _options, _seed, get_patterns(_input, _options)) {}

	bool set_pattern(const Array2D<uint32_t> &pattern, unsigned i, unsigned j) noexcept
	{
		int pattern_id = get_pattern_id(pattern);

		if (pattern_id == -1 || i >= options.get_wave_height() || j >= options.get_wave_width())
		{
			return false;
		}

		set_pattern(pattern_id, i, j);
		return true;
	}

	Array2D<uint32_t> *run() noexcept
	{
		Array2D<unsigned> *result = wfc.run();
		return result == nullptr ? nullptr : to_image(result);
	}
};

enum class Symmetry
{
	X,
	T,
	I,
	L,
	BACKSLASH,
	F
};

constexpr unsigned nb_of_possible_orientations(const Symmetry &symmetry)
{
	return symmetry == Symmetry::X ? 1 : symmetry == Symmetry::I || symmetry == Symmetry::BACKSLASH ? 2
									 : symmetry == Symmetry::T || symmetry == Symmetry::L			? 4
																									: 8;
};

template <typename T>
struct Tile
{
	std::vector<Array2D<T>> data;
	Symmetry symmetry;
	double weight;

	static std::vector<unsigned>
	generate_rotation_map(const Symmetry &symmetry) noexcept
	{
		switch (symmetry)
		{
		case Symmetry::X:
			return {0};
		case Symmetry::I:
		case Symmetry::BACKSLASH:
			return {1, 0};
		case Symmetry::T:
		case Symmetry::L:
			return {1, 2, 3, 0};
		case Symmetry::F:
		default:
			return {1, 2, 3, 0, 5, 6, 7, 4};
		}
	}

	static std::vector<unsigned>
	generate_reflection_map(const Symmetry &symmetry) noexcept
	{
		switch (symmetry)
		{
		case Symmetry::X:
			return {0};
		case Symmetry::I:
			return {0, 1};
		case Symmetry::BACKSLASH:
			return {1, 0};
		case Symmetry::T:
			return {0, 3, 2, 1};
		case Symmetry::L:
			return {1, 0, 3, 2};
		case Symmetry::F:
		default:
			return {4, 7, 6, 5, 0, 3, 2, 1};
		}
	}

	static std::vector<std::vector<unsigned>>
	generate_action_map(const Symmetry &symmetry) noexcept
	{
		std::vector<unsigned> rotation_map = generate_rotation_map(symmetry);
		std::vector<unsigned> reflection_map = generate_reflection_map(symmetry);
		size_t size = rotation_map.size();
		std::vector<std::vector<unsigned>> action_map(8,
													  std::vector<unsigned>(size));
		for (size_t i = 0; i < size; ++i)
		{
			action_map[0][i] = i;
		}

		for (size_t a = 1; a < 4; ++a)
		{
			for (size_t i = 0; i < size; ++i)
			{
				action_map[a][i] = rotation_map[action_map[a - 1][i]];
			}
		}
		for (size_t i = 0; i < size; ++i)
		{
			action_map[4][i] = reflection_map[action_map[0][i]];
		}
		for (size_t a = 5; a < 8; ++a)
		{
			for (size_t i = 0; i < size; ++i)
			{
				action_map[a][i] = rotation_map[action_map[a - 1][i]];
			}
		}
		return action_map;
	}

	static std::vector<Array2D<T>> generate_oriented(Array2D<T> data,
													 Symmetry symmetry) noexcept
	{
		std::vector<Array2D<T>> oriented;
		oriented.push_back(data);

		switch (symmetry)
		{
		case Symmetry::I:
		case Symmetry::BACKSLASH:
			oriented.push_back(data.rotated());
			break;
		case Symmetry::T:
		case Symmetry::L:
			oriented.push_back(data = data.rotated());
			oriented.push_back(data = data.rotated());
			oriented.push_back(data = data.rotated());
			break;
		case Symmetry::F:
			oriented.push_back(data = data.rotated());
			oriented.push_back(data = data.rotated());
			oriented.push_back(data = data.rotated());
			oriented.push_back(data = data.rotated().reflected());
			oriented.push_back(data = data.rotated());
			oriented.push_back(data = data.rotated());
			oriented.push_back(data = data.rotated());
			break;
		default:
			break;
		}

		return oriented;
	}

	Tile(std::vector<Array2D<T>> data, Symmetry symmetry, double weight) noexcept
		: data(data), symmetry(symmetry), weight(weight) {}

	Tile(Array2D<T> data, Symmetry symmetry, double weight) noexcept
		: data(generate_oriented(data, symmetry)), symmetry(symmetry),
		  weight(weight) {}
};

struct TilingWFCOptions
{
	bool periodic_output;
};

template <typename T>
class TilingWFC
{
private:
	std::vector<Tile<T>> tiles;

	std::vector<std::pair<unsigned, unsigned>> id_to_oriented_tile;

	std::vector<std::vector<unsigned>> oriented_tile_ids;

	TilingWFCOptions options;

	WFC wfc;

public:
	unsigned height;
	unsigned width;

private:
	static std::pair<std::vector<std::pair<unsigned, unsigned>>,
					 std::vector<std::vector<unsigned>>>
	generate_oriented_tile_ids(const std::vector<Tile<T>> &tiles) noexcept
	{
		std::vector<std::pair<unsigned, unsigned>> id_to_oriented_tile;
		std::vector<std::vector<unsigned>> oriented_tile_ids;

		unsigned id = 0;
		for (unsigned i = 0; i < tiles.size(); i++)
		{
			oriented_tile_ids.push_back({});
			for (unsigned j = 0; j < tiles[i].data.size(); j++)
			{
				id_to_oriented_tile.push_back({i, j});
				oriented_tile_ids[i].push_back(id);
				id++;
			}
		}

		return {id_to_oriented_tile, oriented_tile_ids};
	}

	static std::vector<std::array<std::vector<unsigned>, 4>> generate_propagator(
		const std::vector<std::tuple<unsigned, unsigned, unsigned, unsigned>>
			&neighbors,
		std::vector<Tile<T>> tiles,
		std::vector<std::pair<unsigned, unsigned>> id_to_oriented_tile,
		std::vector<std::vector<unsigned>> oriented_tile_ids)
	{
		size_t nb_oriented_tiles = id_to_oriented_tile.size();
		std::vector<std::array<std::vector<bool>, 4>> dense_propagator(
			nb_oriented_tiles, {std::vector<bool>(nb_oriented_tiles, false),
								std::vector<bool>(nb_oriented_tiles, false),
								std::vector<bool>(nb_oriented_tiles, false),
								std::vector<bool>(nb_oriented_tiles, false)});

		for (auto neighbor : neighbors)
		{
			unsigned tile1 = std::get<0>(neighbor);
			unsigned orientation1 = std::get<1>(neighbor);
			unsigned tile2 = std::get<2>(neighbor);
			unsigned orientation2 = std::get<3>(neighbor);
			std::vector<std::vector<unsigned>> action_map1 =
				Tile<T>::generate_action_map(tiles[tile1].symmetry);
			std::vector<std::vector<unsigned>> action_map2 =
				Tile<T>::generate_action_map(tiles[tile2].symmetry);

			auto add = [&](unsigned action, unsigned direction)
			{
				unsigned temp_orientation1 = action_map1[action][orientation1];
				unsigned temp_orientation2 = action_map2[action][orientation2];
				unsigned oriented_tile_id1 =
					oriented_tile_ids[tile1][temp_orientation1];
				unsigned oriented_tile_id2 =
					oriented_tile_ids[tile2][temp_orientation2];
				dense_propagator[oriented_tile_id1][direction][oriented_tile_id2] =
					true;
				direction = get_opposite_direction(direction);
				dense_propagator[oriented_tile_id2][direction][oriented_tile_id1] =
					true;
			};

			add(0, 2);
			add(1, 0);
			add(2, 1);
			add(3, 3);
			add(4, 1);
			add(5, 3);
			add(6, 2);
			add(7, 0);
		}

		std::vector<std::array<std::vector<unsigned>, 4>> propagator(
			nb_oriented_tiles);
		for (size_t i = 0; i < nb_oriented_tiles; ++i)
		{
			for (size_t j = 0; j < nb_oriented_tiles; ++j)
			{
				for (size_t d = 0; d < 4; ++d)
				{
					if (dense_propagator[i][d][j])
					{
						propagator[i][d].push_back(j);
					}
				}
			}
		}

		return propagator;
	}

	static std::vector<double>
	get_tiles_weights(const std::vector<Tile<T>> &tiles)
	{
		std::vector<double> frequencies;
		for (size_t i = 0; i < tiles.size(); ++i)
		{
			for (size_t j = 0; j < tiles[i].data.size(); ++j)
			{
				frequencies.push_back(tiles[i].weight / tiles[i].data.size());
			}
		}
		return frequencies;
	}

	Array2D<T> *id_to_tiling(Array2D<unsigned> *ids)
	{
		unsigned size = tiles[0].data[0].height;
		Array2D<T> *tiling = new Array2D<T>(size * ids->height, size * ids->width);
		for (unsigned i = 0; i < ids->height; i++)
		{
			for (unsigned j = 0; j < ids->width; j++)
			{
				std::pair<unsigned, unsigned> oriented_tile =
					id_to_oriented_tile[ids->get(i, j)];
				for (unsigned y = 0; y < size; y++)
				{
					for (unsigned x = 0; x < size; x++)
					{
						tiling->get(i * size + y, j * size + x) =
							tiles[oriented_tile.first].data[oriented_tile.second].get(y, x);
					}
				}
			}
		}
		return tiling;
	}

	void set_tile(unsigned tile_id, unsigned i, unsigned j) noexcept
	{
		for (unsigned p = 0; p < id_to_oriented_tile.size(); p++)
		{
			if (tile_id != p)
			{
				wfc.remove_wave_pattern(i, j, p);
			}
		}
	}

public:
	TilingWFC(
		const std::vector<Tile<T>> &tiles,
		const std::vector<std::tuple<unsigned, unsigned, unsigned, unsigned>>
			&neighbors,
		const unsigned height, const unsigned width,
		const TilingWFCOptions &options, unsigned seed)
		: tiles(tiles),
		  id_to_oriented_tile(generate_oriented_tile_ids(tiles).first),
		  oriented_tile_ids(generate_oriented_tile_ids(tiles).second),
		  options(options),
		  wfc(options.periodic_output, seed, get_tiles_weights(tiles),
			  generate_propagator(neighbors, tiles, id_to_oriented_tile,
								  oriented_tile_ids),
			  height, width),
		  height(height), width(width) {}

	bool set_tile(unsigned tile_id, unsigned orientation, unsigned i, unsigned j) noexcept
	{
		if (tile_id >= oriented_tile_ids.size() || orientation >= oriented_tile_ids[tile_id].size() || i >= height || j >= width)
		{
			return false;
		}

		unsigned oriented_tile_id = oriented_tile_ids[tile_id][orientation];
		set_tile(oriented_tile_id, i, j);
		return true;
	}

	Array2D<T> *run()
	{
		Array2D<unsigned> *a = wfc.run();
		if (a == nullptr)
		{
			return nullptr;
		}
		return id_to_tiling(a);
	}
};

#endif