#pragma once

#include <string>
#include <vector>

namespace parse {

struct Alignment {
    std::vector<std::string> names;
    std::vector<std::string> sequences;
    size_t sites = 0;
};

struct ModelConfig {
    int states = 4;
    std::string subst_model;
    int ncat = 1;
    double alpha = 1.0;
    double pinv = 0.0;
    std::vector<double> freqs;
    std::vector<double> rates;
    std::vector<double> rate_weights;
    bool per_rate_scaling = false;
};

struct FilesConfig {
    std::string alignment;
    std::string tree;
};

struct RunConfig {
    ModelConfig model;
    FilesConfig files;
};

RunConfig parse_config_file(const std::string& path);
Alignment read_alignment_file(const std::string& path);
struct RunInputs {
    RunConfig config;
    Alignment alignment;
    std::string tree;
};

RunInputs load_inputs_from_config(const std::string& config_path);

} // namespace parse
