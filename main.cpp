#include "src/include/overlapping_wfc.h"
#include "src/include/tiling_wfc.h"
#include "src/include/utils/color.h"
#include "src/include/utils/array3D.h"
#include "src/include/wfc.h"
#include "lib/all.h"

#include <unordered_set>
#include <iostream>
#include <chrono>

/**
 * Get a random seed.
 * This function use random_device on linux, but use the C rand function for
 * other targets. This is because, for instance, windows don't implement
 * random_device non-deterministically.
 */
int get_random_seed()
{
#ifdef __linux__
	return random_device()();
#else
	return rand();
#endif
}

Array2D<Color> *read_image(const std::string &file_path) noexcept
{
	int width;
	int height;
	int num_components;
	unsigned char *data = stbi_load(file_path.c_str(), &width, &height, &num_components, 3);
	if (data == nullptr)
	{
		return nullptr;
	}
	Array2D<Color> *m = new Array2D<Color>(height, width);
	for (unsigned i = 0; i < (unsigned)height; i++)
	{
		for (unsigned j = 0; j < (unsigned)width; j++)
		{
			unsigned index = 3 * (i * width + j);
			m->data[i * width + j] = {data[index], data[index + 1], data[index + 2]};
		}
	}
	free(data);
	return m;
}

/**
 * Write an image in the png format.
 */
void write_image_png(const std::string &file_path, const Array2D<Color> &m) noexcept
{
	stbi_write_png(file_path.c_str(), m.width, m.height, 3, (const unsigned char *)m.data.data(), 0);
}

/**
 * Read the overlapping wfc problem from the xml node.
 */
void read_overlapping_instance(pugi::xml_node &node)
{
	std::string name = node.attribute("name").as_string();
	unsigned N = node.attribute("N").as_uint();
	bool periodic_output = node.attribute("periodic").as_bool();
	bool periodic_input = node.attribute("periodicInput").as_bool();
	bool ground = node.attribute("ground").as_bool();
	unsigned symmetry = node.attribute("symmetry").as_uint() != 0 ? node.attribute("symmetry").as_uint() : 8;
	unsigned screenshots = 1;
	unsigned width = node.attribute("width").as_uint() != 0 ? node.attribute("width").as_uint() : 48;
	unsigned height = node.attribute("height").as_uint() != 0 ? node.attribute("height").as_uint() : 48;

	std::cout << name << " started!" << std::endl;
	// Stop hardcoding samples
	const std::string image_path = "inputs/" + name + ".png";
	Array2D<Color> *m = read_image(image_path);
	if (m == nullptr)
	{
		throw "Error while loading " + image_path;
	}
	OverlappingWFCOptions options = {
		periodic_input, periodic_output, height, width, symmetry, ground, N};
	for (unsigned i = 0; i < screenshots; i++)
	{
		for (unsigned test = 0; test < 10; test++)
		{
			int seed = get_random_seed();
			OverlappingWFC<Color> wfc(*m, options, seed);
			Array2D<Color> *success = wfc.run();
			if (success != nullptr)
			{
				write_image_png("outputs/" + name + std::to_string(i) + ".png", *success);
				std::cout << name << " finished!" << std::endl;
				break;
			}
			else
			{
				std::cout << "failed!" << std::endl;
			}
		}
	}
}

/**
 * Transform a symmetry name into its Symmetry enum
 */
Symmetry to_symmetry(const std::string &symmetry_name)
{
	if (symmetry_name == "X")
	{
		return Symmetry::X;
	}
	if (symmetry_name == "T")
	{
		return Symmetry::T;
	}
	if (symmetry_name == "I")
	{
		return Symmetry::I;
	}
	if (symmetry_name == "L")
	{
		return Symmetry::L;
	}
	if (symmetry_name == "\\")
	{
		return Symmetry::backslash;
	}
	if (symmetry_name == "F")
	{
		return Symmetry::F;
	}
	return Symmetry::X;
}

/**
 * Read the names of the tiles in the subset in a tiling WFC problem
 */
std::unordered_set<std::string> read_subset_names(pugi::xml_node &node, const std::string &_subset)
{
	std::unordered_set<std::string> subset_names;
	pugi::xml_node subsets_node = node.child("subsets");
	for (pugi::xml_node &subset : subsets_node.children("subset"))
	{
		if (subset.attribute("name").as_string() == _subset)
		{
			for (pugi::xml_node &tile : subsets_node.children("tile"))
			{
				subset_names.insert(tile.attribute("name").as_string());
			}
		}
	}
	return subset_names;
}

/**
 * Read all tiles for a tiling problem
 */
std::unordered_map<std::string, Tile<Color>> read_tiles(pugi::xml_node &root_node, const std::string &current_dir, const std::string &subset, unsigned size)
{
	std::unordered_set<std::string> subset_names = read_subset_names(root_node, subset);
	std::unordered_map<std::string, Tile<Color>> tiles;
	pugi::xml_node tiles_node = root_node.child("tiles");
	for (pugi::xml_node &tile : tiles_node.children("tile"))
	{
		std::string name = tile.attribute("name").as_string();
		if (subset_names.find(name) == subset_names.end())
			continue;
		Symmetry symmetry = to_symmetry(tile.attribute("symmetry").as_string());
		double weight = tile.attribute("weight").as_double() != 0.0 ? tile.attribute("weight").as_double() : 1.0;
		const std::string image_path = current_dir + "/" + name + ".png";
		Array2D<Color> *image = read_image(image_path);

		if (image == nullptr)
		{
			std::vector<Array2D<Color>> images;
			for (unsigned i = 0; i < nb_of_possible_orientations(symmetry); i++)
			{
				const std::string image_path = current_dir + "/" + name + " " + std::to_string(i) + ".png";
				image = read_image(image_path);
				if (image == nullptr)
				{
					throw "Error while loading " + image_path;
				}
				if ((image->width != size) || (image->height != size))
				{
					throw "Image " + image_path + " has wrond size";
				}
				images.push_back(*image);
			}
			Tile<Color> tile = {images, symmetry, weight};
			tiles.insert({name, tile});
		}
		else
		{
			if ((image->width != size) || (image->height != size))
			{
				throw "Image " + image_path + " has wrong size";
			}

			Tile<Color> tile(*image, symmetry, weight);
			tiles.insert({name, tile});
		}
	}

	return tiles;
}

/**
 * Read the neighbors constraints for a tiling problem.
 * A value {t1,o1,t2,o2} means that the tile t1 with orientation o1 can be
 * placed at the right of the tile t2 with orientation o2.
 */
std::vector<std::tuple<std::string, unsigned, std::string, unsigned>> read_neighbors(pugi::xml_node &root_node)
{
	std::vector<std::tuple<std::string, unsigned, std::string, unsigned>> neighbors;
	pugi::xml_node neighbor_node = root_node.child("neighbors");
	for (pugi::xml_node &neighbor : neighbor_node)
	{
		std::string left = neighbor.attribute("left").as_string();
		std::string::size_type left_delimiter = left.find(" ");
		std::string left_tile = left.substr(0, left_delimiter);
		unsigned left_orientation = 0;
		if (left_delimiter != std::string::npos)
		{
			left_orientation = stoi(left.substr(left_delimiter, std::string::npos));
		}

		std::string right = neighbor.attribute("right").as_string();
		std::string::size_type right_delimiter = right.find(" ");
		std::string right_tile = right.substr(0, right_delimiter);
		unsigned right_orientation = 0;
		if (right_delimiter != std::string::npos)
		{
			right_orientation = stoi(right.substr(right_delimiter, std::string::npos));
		}
		neighbors.push_back(
			{left_tile, left_orientation, right_tile, right_orientation});
	}
	return neighbors;
}

/**
 * Read an instance of a tiling WFC problem.
 */
void read_simpletiled_instance(pugi::xml_node &node, const std::string &current_dir) noexcept
{
	std::string name = node.attribute("name").as_string();
	std::string subset = node.attribute("subset").as_string();
	bool periodic_output = node.attribute("periodic").as_bool();
	unsigned width = node.attribute("width").as_uint() != 0 ? node.attribute("width").as_uint() : 48;
	unsigned height = node.attribute("height").as_uint() != 0 ? node.attribute("height").as_uint() : 48;

	std::cout << name << " " << subset << " started!" << std::endl;

	pugi::xml_document doc;
	std::string path = "inputs/" + name + "/data.xml";
	doc.load_file(path.c_str());
	pugi::xml_node set_node = doc.child("set");
	unsigned size = set_node.attribute("size").as_uint();

	std::unordered_map<std::string, Tile<Color>> tiles_map = read_tiles(set_node, current_dir + "/" + name, subset, size);
	std::unordered_map<std::string, unsigned> tiles_id;
	std::vector<Tile<Color>> tiles;
	unsigned id = 0;
	for (std::pair<std::string, Tile<Color>> tile : tiles_map)
	{
		tiles_id.insert({tile.first, id});
		tiles.push_back(tile.second);
		id++;
	}

	std::vector<std::tuple<std::string, unsigned, std::string, unsigned>> neighbors = read_neighbors(set_node);
	std::vector<std::tuple<unsigned, unsigned, unsigned, unsigned>> neighbors_ids;
	for (auto neighbor : neighbors)
	{
		const std::string &neighbor1 = std::get<0>(neighbor);
		const int &orientation1 = std::get<1>(neighbor);
		const std::string &neighbor2 = std::get<2>(neighbor);
		const int &orientation2 = std::get<3>(neighbor);
		if (tiles_id.find(neighbor1) == tiles_id.end())
		{
			continue;
		}
		if (tiles_id.find(neighbor2) == tiles_id.end())
		{
			continue;
		}
		neighbors_ids.push_back(std::make_tuple(tiles_id[neighbor1], orientation1,
										   tiles_id[neighbor2], orientation2));
	}

	for (unsigned test = 0; test < 10; test++)
	{
		int seed = get_random_seed();
		TilingWFC<Color> wfc(tiles, neighbors_ids, height, width, {periodic_output},
							 seed);

		// For the summer tileset, place water on the borders, and land in the middle
		if (name == "Summer")
		{
			for (int i = 0; i < height; i++)
			{
				wfc.set_tile(tiles_id["water_a"], 0, i, 0);
				wfc.set_tile(tiles_id["water_a"], 0, i, width - 1);
			}
			for (int j = 0; j < width; j++)
			{
				wfc.set_tile(tiles_id["water_a"], 0, 0, j);
				wfc.set_tile(tiles_id["water_a"], 0, height - 1, j);
			}
			wfc.set_tile(tiles_id["grass"], 0, width / 2, height / 2);
		}

		Array2D<Color> *success = wfc.run();
		if (success != nullptr)
		{
			write_image_png("results/" + name + "_" + subset + ".png", *success);
			std::cout << name << " finished!" << std::endl;
			break;
		}
		else
		{
			std::cout << "failed!" << std::endl;
		}
	}
}

std::vector<std::string> get_files_from_directory(const char *path)
{
	tinydir_dir dir;
	std::vector<std::string> res;

	if (tinydir_open(&dir, path) == -1)
	{
		std::cout << "Error opening directory" << std::endl;
		tinydir_close(&dir);
		return res;
	}
	while (dir.has_next)
	{
		tinydir_file file;
		if (tinydir_readfile(&dir, &file) == -1)
		{
			std::cout << "Error reading file" << std::endl;
			tinydir_close(&dir);
		}
		std::string name(file.name);
		if (name != "." || name != "..")
		{
			std::string full_path = std::string(path) + "/" + name;
			res.push_back(full_path);
		}

		if (tinydir_next(&dir) == -1)
		{
			std::cout << "No more files" << std::endl;
			tinydir_close(&dir);
		}
	}
	return res;
}

void remove_from_directory(const char *path)
{
	std::vector<std::string> files = get_files_from_directory(path);
	for (unsigned i = 0; i < files.size(); i++)
	{
		remove(files[i].c_str());
	}
}

/**
 * Read a configuration file containing multiple wfc problems
 */
void read_config_file() noexcept
{
	remove_from_directory("outputs");
	pugi::xml_document doc;
	doc.load_file("samples.xml");
	pugi::xml_node samples = doc.child("samples");
	for (pugi::xml_node node : samples.children())
	{
		if (std::string(node.name()) == "overlapping")
			read_overlapping_instance(node);
		else
			read_simpletiled_instance(node, "inputs");
	}
}

int main()
{

// Initialize rand for non-linux targets
#ifndef __linux__
	srand(time(nullptr));
#endif

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	read_config_file();

	end = std::chrono::system_clock::now();
	int elapsed_s =
		std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
	int elapsed_ms =
		std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
			.count();
	std::cout << "All samples done in " << elapsed_s << "s, " << elapsed_ms % 1000
			  << "ms.\n";
}
