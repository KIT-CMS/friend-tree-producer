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
build_pipeline_mapping(std::vector<std::string> quantities) {
  // function used to build a map with the structure
  // { "pipeline" : [ "quantity", "quantity_1", "quantity_2", ... ] }
  // containing all the quantities available for a given pipeline
  // first, loop over all quantities and split the string by __ to get all
  // pipelines
  std::map<std::string, std::vector<std::string>> map;
  for (auto quantity : quantities) {
    std::vector<std::string> pipeline_split;
    boost::algorithm::split_regex(pipeline_split, quantity, boost::regex("__"));
    // no match
    if (pipeline_split.size() < 2) {
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

std::map<std::string, std::string>
build_quantities_map(std::string pipeline, std::vector<std::string> quantities,
                     std::vector<std::string> required_quantities,
                     bool &run_required, bool verbose = false) {
  // This function builds a map with the structure
  // { "quantity" : "quantitiy_name_to_be_used" }
  // depending on the pipeline, shifted quantities have to be used if available
  std::map<std::string, std::string> quantity_names;
  // if we are not doing nominal, we have to potentially modify variable names
  if (pipeline != "nominal") {
    std::map<std::string, std::vector<std::string>> pipeline_mapping =
        build_pipeline_mapping(quantities);

    // now loop thought the required quantities
    for (size_t i = 0; i < required_quantities.size(); i++) {
      // check if the quantity is in the pipeline mapping map
      if (pipeline_mapping.find(required_quantities[i]) !=
          pipeline_mapping.end()) {
        // now if the pipeline is found in the vector of the mapping, add the
        // quantity to the list of quantity names
        if (std::find(pipeline_mapping[required_quantities[i]].begin(),
                      pipeline_mapping[required_quantities[i]].end(),
                      pipeline) !=
            pipeline_mapping[required_quantities[i]].end()) {
          std::string shifted_name =
              required_quantities[i] + std::string("__") + pipeline;
          run_required = true;
          quantity_names.insert(std::pair<std::string, std::string>(
              required_quantities[i], shifted_name));

          continue;
        } else {
          quantity_names.insert(std::pair<std::string, std::string>(
              required_quantities[i], required_quantities[i]));
        }
      } else {
        quantity_names.insert(std::pair<std::string, std::string>(
            required_quantities[i], required_quantities[i]));
      }
    }
  } else {
    run_required = true;
    for (size_t i = 0; i < required_quantities.size(); i++) {
      quantity_names.insert(std::pair<std::string, std::string>(
          required_quantities[i], required_quantities[i]));
    }
  }
  // if the lengpph of quantites is not equal to the length of
  // required_quantities, something went wrong
  if (quantity_names.size() != required_quantities.size()) {
    std::cout << "Something went wrong, the number of quantities to be used is "
                 "not equal to the number of required quantities"
              << std::endl;
    // print the two vectors
    std::cout << "The quantities to be used are:" << std::endl;
    for (auto quantity : quantity_names) {
      std::cout << quantity.first << std::endl;
    }
    std::cout << "The required quantities are:" << std::endl;
    for (auto quantity : required_quantities) {
      std::cout << quantity << std::endl;
    }
    exit(1);
  }
  if (verbose) {
    // now print the quantity_names map
    std::cout << "Final set of variables to be used" << std::endl;
    for (auto quantity : quantity_names) {
      std::cout << quantity.first << " " << quantity.second << std::endl;
    }
  }
  return quantity_names;
}

std::vector<std::string>
find_quantities(TTree *inputtree,
                std::vector<std::string> required_quantities) {
  // determine the correct quantities to be used
  std::cout << "Checking required quantities" << std::endl;
  std::vector<std::string> quantities;
  // loop trough the leaves
  TObjArray *leavescopy = inputtree->GetListOfLeaves();
  // print all leaves
  for (int i = 0; i < leavescopy->GetEntries(); i++) {
    std::cout << leavescopy->At(i)->GetName() << std::endl;
  }
  int nLeaves = leavescopy->GetEntries();
  for (int i = 0; i < nLeaves; i++) {
    for (auto quantity : required_quantities) {
      std::string leaf_name = leavescopy->At(i)->GetName();
      if (leaf_name.find(quantity) != std::string::npos) {
        quantities.push_back(leaf_name);
        break;
      }
    }
  }
  return quantities;
}

std::vector<std::string>
find_quantities(std::vector<std::string> available_quantities,
                std::vector<std::string> required_quantities) {
  // determine the correct quantities to be used
  std::cout << "Checking required quantities" << std::endl;
  std::vector<std::string> quantities;
  // loop trough the aviailable quantities
  for (auto available_quantity : available_quantities) {
    for (auto quantity : required_quantities) {
      if (available_quantity.find(quantity) != std::string::npos) {
        quantities.push_back(available_quantity);
        break;
      }
    }
  }
  return quantities;
}

const auto default_float = -10.f;
