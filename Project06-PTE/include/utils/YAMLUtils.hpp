#pragma once

#include <yaml-cpp/yaml.h>

#include <vector>

template <typename T>
std::vector<T> get_yaml_values(YAML::Node& node) {
  auto type = node.Type();
  // single value
  if (type == YAML::NodeType::Scalar) {
    vector<T> ret = {node.as<T>()};
    return ret;
  } else if (type == YAML::NodeType::Sequence) {
    return node.as<vector<T>>();
  } else if (type == YAML::NodeType::Map) {
    auto start = node["start"].as<T>();
    auto end = node["end"].as<T>();
    auto N = node["N"].as<size_t>();
    return linspace<T>(start, end, N);
  } else {
    throw std::runtime_error(
        "Range node must either be list, scalar or a map defining "
        "'start','end' and 'N'.");
  }
}