#ifndef FAST_WFC_OVERLAPPING_WFC_HPP
#define FAST_WFC_OVERLAPPING_WFC_HPP

#include "utils/array2D.hpp"
#include "wfc.hpp"

#include <vector>
#include <algorithm>
#include <unordered_map>

/**
 * Options needed to use the overlapping wfc.
 */
struct OverlappingWFCOptions
{
	bool periodic_input;   // True if the input is toric.
	bool periodic_output;  // True if the output is toric.
	unsigned out_height;   // The height of the output in pixels.
	unsigned out_width;	   // The width of the output in pixels.
	unsigned symmetry;	   // The number of symmetries (the order is defined in wfc).
	bool ground;		   // True if the ground needs to be set (see init_ground).
	unsigned pattern_size; // The width and height in pixel of the patterns.

	/**
	 * Get the wave height given these options.
	 */
	unsigned get_wave_height() const noexcept
	{
		return periodic_output ? out_height : out_height - pattern_size + 1;
	}

	/**
	 * Get the wave width given these options.
	 */
	unsigned get_wave_width() const noexcept
	{
		return periodic_output ? out_width : out_width - pattern_size + 1;
	}
};

/**
 * Class generating a new image with the overlapping WFC algorithm.
 */
class OverlappingWFC
{

private:
	/**
	 * The input image. T is usually a uint_32.
	 */
	Array2D<uint32_t> input;

	/**
	 * Options needed by the algorithm.
	 */
	OverlappingWFCOptions options;

	/**
	 * The array of the different patterns extracted from the input.
	 */
	std::vector<Array2D<uint32_t>> patterns;

	/**
	 * The underlying generic WFC algorithm.
	 */
	WFC wfc;

	/**
	 * Constructor initializing the wfc.
	 * This constructor is called by the other constructors.
	 * This is necessary in order to initialize wfc only once.
	 */
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
		// If necessary, the ground is set.
		if (_options.ground)
		{
			init_ground(wfc, _input, _patterns.first, _options);
		}
	}

	/**
	 * Constructor used only to call the other constructor with more computed
	 * parameters.
	 */
	OverlappingWFC(
		const Array2D<uint32_t> &_input,
		const OverlappingWFCOptions &_options,
		const unsigned &_seed,
		const std::pair<std::vector<Array2D<uint32_t>>, std::vector<double>> &_patterns) noexcept
		: OverlappingWFC(
			  _input, _options, _seed, _patterns,
			  generate_compatible(_patterns.first)) {}

	/**
	 * Init the ground of the output image.
	 * The lowest middle pattern is used as a floor (and ceiling when the input is
	 * toric) and is placed at the lowest possible pattern position in the output
	 * image, on all its width. The pattern cannot be used at any other place in
	 * the output image.
	 */
	void init_ground(
		WFC &_wfc,
		const Array2D<uint32_t> &_input,
		const std::vector<Array2D<uint32_t>> &_patterns,
		const OverlappingWFCOptions &_options) noexcept
	{
		unsigned ground_pattern_id = get_ground_pattern_id(_input, _patterns, _options);

		// Place the pattern in the ground.
		for (unsigned j = 0; j < _options.get_wave_width(); j++)
		{
			set_pattern(ground_pattern_id, _options.get_wave_height() - 1, j);
		}

		// Remove the pattern from the other positions.
		for (unsigned i = 0; i < _options.get_wave_height() - 1; i++)
		{
			for (unsigned j = 0; j < _options.get_wave_width(); j++)
			{
				wfc.remove_wave_pattern(i, j, ground_pattern_id);
			}
		}

		// Propagate the information with wfc.
		_wfc.propagate();
	}

	/**
	 * Return the id of the lowest middle pattern.
	 */
	static unsigned
	get_ground_pattern_id(const Array2D<uint32_t> &_input,
						  const std::vector<Array2D<uint32_t>> &_patterns,
						  const OverlappingWFCOptions &_options) noexcept
	{
		// Get the pattern.
		std::size_t y = _input.height - 1;
		std::size_t x = _input.width / 2;
		unsigned size = _options.pattern_size;
		Array2D<uint32_t> ground_pattern = _input.get_sub_array(y, x, size, size);

		// Retrieve the id of the pattern.
		for (unsigned i = 0; i < _patterns.size(); i++)
		{
			if (ground_pattern == _patterns[i])
			{
				return i;
			}
		}

		// The pattern exists.
		assert(false);
		return 0;
	}

	/**
	 * Return the list of patterns, as well as their probabilities of apparition.
	 */
	static std::pair<std::vector<Array2D<uint32_t>>, std::vector<double>>
	get_patterns(
		const Array2D<uint32_t> &_input,
		const OverlappingWFCOptions &_options) noexcept
	{
		std::vector<Array2D<uint32_t>> _patterns;
		std::unordered_map<Array2D<uint32_t>, unsigned> patterns_id;

		// The number of time a pattern is seen in the input image.
		std::vector<double> patterns_weight;

		unsigned size = _options.pattern_size;
		std::vector<Array2D<uint32_t>> symmetries(8, Array2D<uint32_t>(size, size));

		unsigned max_i = _options.periodic_input ? _input.height : _input.height - size + 1;
		unsigned max_j = _options.periodic_input ? _input.width : _input.width - size + 1;

		// patterns.reserve(max_i * max_j * 8);
		// patterns_weight.reserve(max_i * max_j * 8);
		// unsigned patterns_cont = 0;

		for (unsigned i = 0; i < max_i; i++)
		{
			for (unsigned j = 0; j < max_j; j++)
			{
				// Compute the symmetries of every pattern in the image.
				symmetries[0].data = _input.get_sub_array(i, j, size, size).data;
				symmetries[1].data = symmetries[0].reflected().data;
				symmetries[2].data = symmetries[0].rotated().data;
				symmetries[3].data = symmetries[2].reflected().data;
				symmetries[4].data = symmetries[2].rotated().data;
				symmetries[5].data = symmetries[4].reflected().data;
				symmetries[6].data = symmetries[4].rotated().data;
				symmetries[7].data = symmetries[6].reflected().data;

				// The number of symmetries in the option class define which symetries
				// will be used.
				for (unsigned k = 0; k < _options.symmetry; k++)
				{
					auto res = patterns_id.insert(
						std::make_pair(symmetries[k], _patterns.size()));

					// If the pattern already exist, we just have to increase its number
					// of appearance.
					if (!res.second)
					{
						patterns_weight[res.first->second] += 1;
					}
					else
					{
						_patterns.push_back(symmetries[k]);
						patterns_weight.push_back(1);
						// patterns_cont++;
					}
				}
			}
		}
		// patterns.resize(patterns_cont);
		// patterns_weight.resize(patterns_weight);

		return {_patterns, patterns_weight};
	}

	/**
	 * Return true if the pattern1 is compatible with pattern2
	 * when pattern2 is at a distance (dy,dx) from pattern1.
	 */
	static bool agrees(const Array2D<uint32_t> &pattern1, const Array2D<uint32_t> &pattern2, int dy, int dx) noexcept
	{
		unsigned xmin = dx < 0 ? 0 : dx;
		unsigned xmax = dx < 0 ? dx + pattern2.width : pattern1.width;
		unsigned ymin = dy < 0 ? 0 : dy;
		unsigned ymax = dy < 0 ? dy + pattern2.height : pattern1.width;

		// Iterate on every pixel contained in the intersection of the two pattern.
		for (unsigned y = ymin; y < ymax; y++)
		{
			for (unsigned x = xmin; x < xmax; x++)
			{
				// Check if the color is the same in the two patterns in that pixel.
				if (pattern1.get(y, x) != pattern2.get(y - dy, x - dx))
				{
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * Precompute the function agrees(pattern1, pattern2, dy, dx).
	 * If agrees(pattern1, pattern2, dy, dx), then compatible[pattern1][direction]
	 * contains pattern2, where direction is the direction defined by (dy, dx)
	 * (see direction.hpp).
	 */
	static std::vector<std::array<std::vector<unsigned>, 4>>
	generate_compatible(const std::vector<Array2D<uint32_t>> &_patterns) noexcept
	{
		std::vector<std::array<std::vector<unsigned>, 4>> compatible =
			std::vector<std::array<std::vector<unsigned>, 4>>(_patterns.size());

		// Iterate on every dy, dx, pattern1 and pattern2
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

	/**
	 * Transform a 2D array containing the patterns id to a 2D array containing
	 * the pixels.
	 */
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

	/**
	 * Set the pattern at a specific position, given its pattern id
	 * pattern_id needs to be a valid pattern id, and i and j needs to be in the wave range
	 */
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
	/**
	 * The constructor used by the user.
	 */
	OverlappingWFC(
		const Array2D<uint32_t> &_input,
		const OverlappingWFCOptions &_options,
		int _seed) noexcept
		: OverlappingWFC(_input, _options, _seed, get_patterns(_input, _options)) {}

	/**
	 * Set the pattern at a specific position.
	 * Returns false if the given pattern does not exist, or if the
	 * coordinates are not in the wave
	 */
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

	/**
	 * Run the WFC algorithm, and return the result if the algorithm succeeded.
	 */
	Array2D<uint32_t> *run() noexcept
	{
		Array2D<unsigned> *result = wfc.run();
		return result == nullptr ? nullptr : to_image(result);
	}
};

#endif // FAST_WFC_WFC_HPP_
