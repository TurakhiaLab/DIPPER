#include "parse_file.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace parse {

namespace detail {

enum class ConfigFormat {
    YAML,
    JSON
};

enum class JsonType {
    Null,
    Object,
    Array,
    String,
    Number,
    Bool
};

struct JsonValue {
    JsonType type = JsonType::Null;
    std::string string_value;
    double number_value = 0.0;
    bool bool_value = false;
    std::vector<std::shared_ptr<JsonValue>> array;
    std::unordered_map<std::string, std::shared_ptr<JsonValue>> object;
};

inline std::string trim_copy(const std::string& s) {
    size_t head = 0;
    while (head < s.size() && std::isspace(static_cast<unsigned char>(s[head]))) ++head;
    size_t tail = s.size();
    while (tail > head && std::isspace(static_cast<unsigned char>(s[tail - 1]))) --tail;
    return s.substr(head, tail - head);
}

inline std::string strip_quotes(const std::string& s) {
    if (s.size() >= 2) {
        if ((s.front() == '"' && s.back() == '"') ||
            (s.front() == '\'' && s.back() == '\'')) {
            return s.substr(1, s.size() - 2);
        }
    }
    return s;
}

inline std::string strip_yaml_comment(const std::string& raw) {
    bool in_double = false;
    bool in_single = false;
    for (size_t i = 0; i < raw.size(); ++i) {
        char c = raw[i];
        if (c == '"' && !in_single) {
            in_double = !in_double;
        } else if (c == '\'' && !in_double) {
            in_single = !in_single;
        }
        if (c == '#' && !in_double && !in_single) {
            return raw.substr(0, i);
        }
    }
    return raw;
}

inline void skip_json_ws(const std::string& text, size_t& pos) {
    while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
}

JsonValue parse_json_value(const std::string& text, size_t& pos);

std::string parse_json_string(const std::string& text, size_t& pos) {
    if (pos >= text.size() || text[pos] != '"')
        throw std::runtime_error("Expected '\"' while parsing JSON string.");
    ++pos;
    std::string output;
    while (pos < text.size()) {
        char c = text[pos++];
        if (c == '"')
            return output;
        if (c == '\\') {
            if (pos >= text.size())
                throw std::runtime_error("Unterminated escape in JSON string.");
            char esc = text[pos++];
            switch (esc) {
                case '"': output.push_back('"'); break;
                case '\\': output.push_back('\\'); break;
                case '/': output.push_back('/'); break;
                case 'n': output.push_back('\n'); break;
                case 'r': output.push_back('\r'); break;
                case 't': output.push_back('\t'); break;
                case 'b': output.push_back('\b'); break;
                case 'f': output.push_back('\f'); break;
                default:
                    throw std::runtime_error("Unsupported escape sequence in JSON string.");
            }
        } else {
            output.push_back(c);
        }
    }
    throw std::runtime_error("Unterminated JSON string.");
}

double parse_json_number(const std::string& text, size_t& pos) {
    size_t start = pos;
    if (pos < text.size() && (text[pos] == '-' || text[pos] == '+')) ++pos;
    while (pos < text.size() && std::isdigit(static_cast<unsigned char>(text[pos]))) ++pos;
    if (pos < text.size() && text[pos] == '.') {
        ++pos;
        while (pos < text.size() && std::isdigit(static_cast<unsigned char>(text[pos]))) ++pos;
    }
    if (pos < text.size() && (text[pos] == 'e' || text[pos] == 'E')) {
        ++pos;
        if (pos < text.size() && (text[pos] == '-' || text[pos] == '+')) ++pos;
        while (pos < text.size() && std::isdigit(static_cast<unsigned char>(text[pos]))) ++pos;
    }
    if (start == pos)
        throw std::runtime_error("Invalid JSON number.");
    return std::stod(text.substr(start, pos - start));
}

JsonValue parse_json_object(const std::string& text, size_t& pos) {
    JsonValue out;
    out.type = JsonType::Object;
    ++pos;
    skip_json_ws(text, pos);
    if (pos < text.size() && text[pos] == '}') {
        ++pos;
        return out;
    }
    while (pos < text.size()) {
        skip_json_ws(text, pos);
        std::string key = parse_json_string(text, pos);
        skip_json_ws(text, pos);
        if (pos >= text.size() || text[pos] != ':')
            throw std::runtime_error("Missing ':' in JSON object.");
        ++pos;
        skip_json_ws(text, pos);
        JsonValue value = parse_json_value(text, pos);
        out.object.emplace(std::move(key), std::make_shared<JsonValue>(std::move(value)));
        skip_json_ws(text, pos);
        if (pos >= text.size())
            throw std::runtime_error("Unexpected end of JSON while parsing object.");
        char c = text[pos];
        if (c == ',') {
            ++pos;
            continue;
        } else if (c == '}') {
            ++pos;
            break;
        } else {
            throw std::runtime_error("Invalid JSON object separator.");
        }
    }
    return out;
}

JsonValue parse_json_array(const std::string& text, size_t& pos) {
    JsonValue out;
    out.type = JsonType::Array;
    ++pos;
    skip_json_ws(text, pos);
    if (pos < text.size() && text[pos] == ']') {
        ++pos;
        return out;
    }
    while (pos < text.size()) {
        skip_json_ws(text, pos);
        JsonValue value = parse_json_value(text, pos);
        out.array.push_back(std::make_shared<JsonValue>(std::move(value)));
        skip_json_ws(text, pos);
        if (pos >= text.size())
            throw std::runtime_error("Unexpected end of JSON while parsing array.");
        char c = text[pos];
        if (c == ',') {
            ++pos;
            continue;
        } else if (c == ']') {
            ++pos;
            break;
        } else {
            throw std::runtime_error("Invalid JSON array separator.");
        }
    }
    return out;
}

JsonValue parse_json_value(const std::string& text, size_t& pos) {
    skip_json_ws(text, pos);
    if (pos >= text.size())
        throw std::runtime_error("Unexpected end of JSON.");
    char c = text[pos];
    if (c == '{') return parse_json_object(text, pos);
    if (c == '[') return parse_json_array(text, pos);
    if (c == '"') {
        JsonValue out;
        out.type = JsonType::String;
        out.string_value = parse_json_string(text, pos);
        return out;
    }
    if (std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+') {
        JsonValue out;
        out.type = JsonType::Number;
        out.number_value = parse_json_number(text, pos);
        return out;
    }
    if (text.compare(pos, 4, "null") == 0) {
        pos += 4;
        JsonValue out;
        out.type = JsonType::Null;
        return out;
    }
    if (text.compare(pos, 4, "true") == 0) {
        pos += 4;
        JsonValue out;
        out.type = JsonType::Bool;
        out.bool_value = true;
        return out;
    }
    if (text.compare(pos, 5, "false") == 0) {
        pos += 5;
        JsonValue out;
        out.type = JsonType::Bool;
        out.bool_value = false;
        return out;
    }
    throw std::runtime_error("Unsupported JSON value.");
}

int json_require_int(const JsonValue& value) {
    if (value.type == JsonType::Number) {
        return static_cast<int>(value.number_value);
    }
    if (value.type == JsonType::String) {
        return std::stoi(value.string_value);
    }
    throw std::runtime_error("Expected numeric value in JSON.");
}

double json_require_double(const JsonValue& value) {
    if (value.type == JsonType::Number) {
        return value.number_value;
    }
    if (value.type == JsonType::String) {
        return std::stod(value.string_value);
    }
    throw std::runtime_error("Expected double value in JSON.");
}

bool json_require_bool(const JsonValue& value) {
    if (value.type == JsonType::Bool)
        return value.bool_value;
    if (value.type == JsonType::String) {
        std::string lower = value.string_value;
        for (char& ch : lower) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        if (lower == "true" || lower == "1") return true;
        if (lower == "false" || lower == "0") return false;
    }
    throw std::runtime_error("Expected boolean value in JSON.");
}

std::string json_require_string(const JsonValue& value) {
    if (value.type == JsonType::String) return value.string_value;
    if (value.type == JsonType::Number) {
        std::ostringstream oss;
        oss << value.number_value;
        return oss.str();
    }
    throw std::runtime_error("Expected string value in JSON.");
}

std::vector<double> json_require_number_array(const JsonValue& value) {
    if (value.type != JsonType::Array)
        throw std::runtime_error("Expected JSON array.");
    std::vector<double> result;
    result.reserve(value.array.size());
    for (const auto& elt : value.array) {
        if (!elt)
            throw std::runtime_error("Null JSON array element.");
        result.push_back(json_require_double(*elt));
    }
    return result;
}

ModelConfig parse_model_json(const JsonValue& obj) {
    if (obj.type != JsonType::Object)
        throw std::runtime_error("JSON 'model' section must be an object.");
    ModelConfig out;
    auto find = obj.object.find("states");
    if (find == obj.object.end()) throw std::runtime_error("Missing 'states' in JSON model.");
    out.states = json_require_int(*find->second);
    auto subst = obj.object.find("subst_model");
    if (subst == obj.object.end()) throw std::runtime_error("Missing 'subst_model' in JSON model.");
    out.subst_model = json_require_string(*subst->second);
    auto ncat = obj.object.find("ncat");
    if (ncat == obj.object.end()) throw std::runtime_error("Missing 'ncat' in JSON model.");
    out.ncat = json_require_int(*ncat->second);
    auto alpha = obj.object.find("alpha");
    if (alpha != obj.object.end()) out.alpha = json_require_double(*alpha->second);
    auto pinv = obj.object.find("pinv");
    if (pinv != obj.object.end()) out.pinv = json_require_double(*pinv->second);
    auto freqs = obj.object.find("freqs");
    if (freqs == obj.object.end()) throw std::runtime_error("Missing 'freqs' in JSON model.");
    out.freqs = json_require_number_array(*freqs->second);
    auto rates = obj.object.find("rates");
    if (rates == obj.object.end()) throw std::runtime_error("Missing 'rates' in JSON model.");
    out.rates = json_require_number_array(*rates->second);
    auto rate_weights = obj.object.find("rate_weights");
    if (rate_weights != obj.object.end()) out.rate_weights = json_require_number_array(*rate_weights->second);
    auto per_rate_scaling = obj.object.find("per_rate_scaling");
    if (per_rate_scaling != obj.object.end())
        out.per_rate_scaling = json_require_bool(*per_rate_scaling->second);
    return out;
}

FilesConfig parse_files_json(const JsonValue& obj) {
    if (obj.type != JsonType::Object)
        throw std::runtime_error("JSON 'files' section must be an object.");
    FilesConfig out;
    auto aln = obj.object.find("alignment");
    if (aln == obj.object.end())
        throw std::runtime_error("Missing 'alignment' in JSON files.");
    out.alignment = json_require_string(*aln->second);
    auto tree = obj.object.find("tree");
    if (tree == obj.object.end())
        throw std::runtime_error("Missing 'tree' in JSON files.");
    out.tree = json_require_string(*tree->second);
    return out;
}

std::optional<std::string> yaml_get_optional_value(
    const std::unordered_map<std::string, std::string>& map,
    const std::string& key)
{
    auto it = map.find(key);
    if (it == map.end())
        return std::nullopt;
    return strip_quotes(trim_copy(it->second));
}

std::string yaml_get_required_value(
    const std::unordered_map<std::string, std::string>& map,
    const std::string& key)
{
    auto it = map.find(key);
    if (it == map.end())
        throw std::runtime_error("Missing '" + key + "' in YAML.");
    return strip_quotes(trim_copy(it->second));
}

int yaml_require_int(
    const std::unordered_map<std::string, std::string>& map,
    const std::string& key)
{
    return std::stoi(yaml_get_required_value(map, key));
}

double yaml_require_double(
    const std::unordered_map<std::string, std::string>& map,
    const std::string& key)
{
    return std::stod(yaml_get_required_value(map, key));
}

bool yaml_parse_bool(const std::string& raw) {
    std::string lower = raw;
    for (char& ch : lower) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (lower == "true" || lower == "1") return true;
    if (lower == "false" || lower == "0") return false;
    throw std::runtime_error("Invalid YAML boolean: " + raw);
}

bool yaml_get_bool(
    const std::unordered_map<std::string, std::string>& map,
    const std::string& key,
    bool default_value)
{
    auto opt = yaml_get_optional_value(map, key);
    if (!opt) return default_value;
    return yaml_parse_bool(*opt);
}

std::vector<double> yaml_parse_array(const std::string& raw) {
    std::string trimmed = trim_copy(raw);
    if (trimmed.empty())
        return {};
    if (trimmed.front() == '[') trimmed.erase(trimmed.begin());
    if (!trimmed.empty() && trimmed.back() == ']') trimmed.pop_back();
    std::vector<double> values;
    std::istringstream iss(trimmed);
    while (iss.good()) {
        std::string token;
        if (!std::getline(iss, token, ',')) break;
        auto item = trim_copy(token);
        if (item.empty()) continue;
        values.push_back(std::stod(item));
    }
    return values;
}

std::unordered_map<std::string, std::string> parse_yaml_section(
    const std::string& text,
    const std::string& section_name)
{
    std::unordered_map<std::string, std::string> kv;
    std::istringstream iss(text);
    std::string line;
    bool in_section = false;
    int base_indent = 0;
    while (std::getline(iss, line)) {
        size_t idx = 0;
        while (idx < line.size() && (line[idx] == ' ' || line[idx] == '\t')) ++idx;
        std::string trimmed = trim_copy(line.substr(idx));
        if (trimmed.empty() || trimmed[0] == '#') continue;
        if (!in_section) {
            size_t colon = trimmed.find(':');
            if (colon == std::string::npos) continue;
            std::string key = trim_copy(trimmed.substr(0, colon));
            if (key == section_name) {
                in_section = true;
                base_indent = static_cast<int>(idx);
                continue;
            }
        } else {
            if (static_cast<int>(idx) <= base_indent) {
                break;
            }
            size_t colon = trimmed.find(':');
            if (colon == std::string::npos) continue;
            std::string key = trim_copy(trimmed.substr(0, colon));
            std::string value = strip_yaml_comment(trimmed.substr(colon + 1));
            kv.emplace(key, value);
        }
    }
    if (!in_section)
        throw std::runtime_error("Missing '" + section_name + "' section in YAML config.");
    return kv;
}

ModelConfig build_model_from_yaml_map(
    const std::unordered_map<std::string, std::string>& map)
{
    ModelConfig out;
    out.states = yaml_require_int(map, "states");
    out.subst_model = yaml_get_required_value(map, "subst_model");
    out.ncat = yaml_require_int(map, "ncat");
    if (auto opt_alpha = yaml_get_optional_value(map, "alpha"))
        out.alpha = std::stod(*opt_alpha);
    if (auto opt_pinv = yaml_get_optional_value(map, "pinv"))
        out.pinv = std::stod(*opt_pinv);
    out.freqs = yaml_parse_array(yaml_get_required_value(map, "freqs"));
    out.rates = yaml_parse_array(yaml_get_required_value(map, "rates"));
    auto maybe_weights = yaml_get_optional_value(map, "rate_weights");
    if (maybe_weights) {
        out.rate_weights = yaml_parse_array(*maybe_weights);
    }
    out.per_rate_scaling = yaml_get_bool(map, "per_rate_scaling", false);
    return out;
}

FilesConfig build_files_from_yaml_map(
    const std::unordered_map<std::string, std::string>& map)
{
    FilesConfig out;
    out.alignment = yaml_get_required_value(map, "alignment");
    out.tree = yaml_get_required_value(map, "tree");
    return out;
}

RunConfig parse_yaml_config(const std::string& text) {
    auto model_map = parse_yaml_section(text, "model");
    auto file_map = parse_yaml_section(text, "files");
    return RunConfig{ build_model_from_yaml_map(model_map), build_files_from_yaml_map(file_map) };
}

RunConfig parse_json_config(const std::string& text) {
    size_t pos = 0;
    JsonValue root = parse_json_value(text, pos);
    if (root.type != JsonType::Object)
        throw std::runtime_error("Top-level JSON config must be an object.");
    auto model_it = root.object.find("model");
    if (model_it == root.object.end())
        throw std::runtime_error("Missing 'model' section in JSON config.");
    auto files_it = root.object.find("files");
    if (files_it == root.object.end())
        throw std::runtime_error("Missing 'files' section in JSON config.");
    return RunConfig{ parse_model_json(*model_it->second), parse_files_json(*files_it->second) };
}

ConfigFormat detect_config_format(const std::string& text) {
    size_t idx = 0;
    while (idx < text.size() && std::isspace(static_cast<unsigned char>(text[idx]))) ++idx;
    if (idx < text.size() && text[idx] == '{')
        return ConfigFormat::JSON;
    return ConfigFormat::YAML;
}

Alignment read_fasta_alignment(const std::string& path) {
    Alignment aln;
    std::ifstream ifs(path);
    if (!ifs)
        throw std::runtime_error("Cannot open alignment file: " + path);
    std::string line;
    std::string current_name;
    std::string current_seq;
    auto push_current = [&]() {
        if (!current_name.empty()) {
            aln.names.push_back(current_name);
            aln.sequences.push_back(current_seq);
            current_seq.clear();
        }
    };
    while (std::getline(ifs, line)) {
        std::string trimmed = trim_copy(line);
        if (trimmed.empty()) continue;
        if (trimmed[0] == '>') {
            push_current();
            current_name = trim_copy(trimmed.substr(1));
        } else if (!current_name.empty()) {
            for (char c : trimmed) {
                if (!std::isspace(static_cast<unsigned char>(c)))
                    current_seq.push_back(c);
            }
        }
    }
    push_current();
    if (aln.names.empty())
        throw std::runtime_error("FASTA alignment has no sequences.");
    aln.sites = aln.sequences[0].size();
    for (const auto& seq : aln.sequences) {
        if (seq.size() != aln.sites)
            throw std::runtime_error("Inconsistent sequence lengths in alignment.");
    }
    return aln;
}

std::vector<std::string> split_on_whitespace(const std::string& text) {
    std::vector<std::string> tokens;
    std::istringstream iss(text);
    std::string token;
    while (iss >> token) tokens.push_back(token);
    return tokens;
}

Alignment read_phylip_alignment(const std::string& path) {
    Alignment aln;
    std::ifstream ifs(path);
    if (!ifs)
        throw std::runtime_error("Cannot open alignment file: " + path);
    std::string header;
    while (std::getline(ifs, header)) {
        if (trim_copy(header).empty()) continue;
        break;
    }
    if (header.empty())
        throw std::runtime_error("PHYLIP file is empty.");
    auto tokens = split_on_whitespace(trim_copy(header));
    if (tokens.size() < 2)
        throw std::runtime_error("PHYLIP header must contain sequence count and site count.");
    int sequence_count = std::stoi(tokens[0]);
    int site_count = std::stoi(tokens[1]);
    std::string line;
    while (aln.names.size() < static_cast<size_t>(sequence_count) && std::getline(ifs, line)) {
        std::string trimmed = trim_copy(line);
        if (trimmed.empty()) continue;
        auto parts = split_on_whitespace(trimmed);
        if (parts.empty()) continue;
        std::string name = parts[0];
        std::string seq;
        for (size_t i = 1; i < parts.size(); ++i)
            seq += parts[i];
        if (seq.empty())
            throw std::runtime_error("Empty sequence encountered in PHYLIP.");
        if (seq.size() != static_cast<size_t>(site_count))
            throw std::runtime_error("PHYLIP sequence length mismatch header.");
        aln.names.push_back(name);
        aln.sequences.push_back(seq);
    }
    if (aln.names.size() != static_cast<size_t>(sequence_count))
        throw std::runtime_error("PHYLIP sequence count mismatch.");
    aln.sites = site_count;
    return aln;
}

Alignment read_alignment(const std::string& path) {
    std::ifstream probe(path);
    if (!probe)
        throw std::runtime_error("Cannot open alignment file: " + path);
    std::string first_line;
    if (!std::getline(probe, first_line))
        throw std::runtime_error("Alignment file is empty.");
    std::string trimmed = trim_copy(first_line);
    probe.close();
    if (trimmed.empty())
        throw std::runtime_error("Alignment file is empty.");
    if (trimmed.front() == '>')
        return read_fasta_alignment(path);
    return read_phylip_alignment(path);
}

} // namespace detail

RunConfig parse_config_file(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs)
        throw std::runtime_error("Cannot open config file: " + path);
    std::string all;
    {
        std::ostringstream oss;
        oss << ifs.rdbuf();
        all = oss.str();
    }
    detail::ConfigFormat format = detail::detect_config_format(all);
    if (format == detail::ConfigFormat::JSON)
        return detail::parse_json_config(all);
    return detail::parse_yaml_config(all);
}

Alignment read_alignment_file(const std::string& path) {
    return detail::read_alignment(path);
}

RunInputs load_inputs_from_config(const std::string& config_path) {
    RunConfig config = parse_config_file(config_path);

    namespace fs = std::filesystem;
    auto make_absolute = [&](const std::string& target) {
        fs::path candidate(target);
        if (candidate.is_absolute()) return candidate;
        fs::path base = fs::path(config_path).parent_path();
        if (base.empty()) base = fs::current_path();
        return base / candidate;
    };

    Alignment alignment = detail::read_alignment(make_absolute(config.files.alignment).string());

    std::ifstream tree_stream(make_absolute(config.files.tree).string());
    if (!tree_stream)
        throw std::runtime_error("Cannot open tree file: " + config.files.tree);
    std::ostringstream tree_buf;
    tree_buf << tree_stream.rdbuf();
    std::string tree_text = tree_buf.str();

    return RunInputs{ std::move(config), std::move(alignment), std::move(tree_text) };
}

} // namespace parse
