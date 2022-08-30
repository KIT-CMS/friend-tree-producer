#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string_regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

namespace fs = boost::filesystem;

std::string folder_to_channel(std::string foldername) {
  std::vector<std::string> foldername_split;
  boost::split(foldername_split, foldername, [](char c) { return c == '_'; });
  return foldername_split.at(0);
}

std::string filename_from_inputpath(std::string input) {
  std::vector<std::string> path_split;
  boost::split(path_split, input, [](char c) { return c == '/'; });
  std::string filename = path_split.at(path_split.size() - 1);
  boost::replace_all(filename, ".root", "");
  return filename;
}

std::string outputname_from_settings(std::string input, std::string folder,
                                     unsigned int first_entry,
                                     unsigned int last_entry,
                                     fs::path output_dir = "",
                                     bool organize_outputs = true) {
  std::string filename = filename_from_inputpath(input);
  fs::path full_path = "";
  if (organize_outputs) {
    std::string outputname = filename + "/" + filename + "_" + folder + "_" +
                             std::to_string(first_entry) + "_" +
                             std::to_string(last_entry) + ".root";
    full_path = output_dir / outputname;
  } else {
    std::string outputname = filename + "_" + folder + "_" +
                             std::to_string(first_entry) + "_" +
                             std::to_string(last_entry) + ".root";
    full_path = outputname;
  }
  return full_path.string();
}

std::string outputname_from_settings_crown(
    std::string input, std::string folder, unsigned int first_entry,
    unsigned int last_entry, fs::path output_dir = "", std::string era = "2018",
    std::string channel = "tt", bool organize_outputs = true) {
  std::cout << "Using input " << input << " folder " << folder
            << " first entry " << first_entry << " last entry " << last_entry
            << " era " << era << " channel " << channel << std::endl;
  std::string filename = filename_from_inputpath(input);
  fs::path full_path = "";
  if (organize_outputs) {
    std::string outputname = era + "_" + channel + "_" + filename + "/" +
                             filename + "_" + folder + "_" + era + "_" +
                             channel + "_" + std::to_string(first_entry) + "_" +
                             std::to_string(last_entry) + ".root";
    full_path = output_dir / outputname;
  } else {
    std::string outputname = filename + "_" + folder + "_" + era + "_" +
                             channel + "_" + std::to_string(first_entry) + "_" +
                             std::to_string(last_entry) + ".root";
    full_path = outputname;
  }
  return full_path.string();
}

std::map<std::string, std::vector<std::string>>
build_quantities_map(std::vector<std::string> quantities) {
  // first, loop over all quantities and split the string by __ to get all
  // pipelines
  std::map<std::string, std::vector<std::string>> map;
  for (auto quantity : quantities) {
    std::vector<std::string> pipeline_split;
    boost::algorithm::split_regex(pipeline_split, quantity, boost::regex("__"));
    // boost::split(pipeline_split,quantity, boost::is_any_of("__"));
    if (pipeline_split.size() < 2) {
      std::cout << "Warning: quantity " << quantity
                << " does not contain a pipeline, cont." << std::endl;
      continue;
    }
    std::string quantity_name = pipeline_split.at(0);
    std::string pipeline = pipeline_split.at(1);
    if (map.find(quantity_name) == map.end()) {
      map[quantity_name] = std::vector<std::string>();
    }
    map[quantity_name].push_back(pipeline);
  }
  return map;
}

const auto default_float = -10.f;
