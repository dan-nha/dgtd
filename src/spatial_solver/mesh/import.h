#ifndef IMPORT_H
#define IMPORT_H

#include "mesh.h"

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Import {

/**
 * Import n-th entry of a line
 *
 * @param[in] line Line to read from
 * @param[in] entry_number entry number in line
 *
 * @return The n-th entry of a line
 */
template <typename T>
const T get_entry(const std::string line, const size_t entry_number = 1) {

  if (entry_number == 0)
    return T(0);

  std::stringstream sstr(line);
  T output;

  for (size_t n = 0; n < entry_number; ++n)
    sstr >> output;

  return output;
}

/**
 * Import a number of entries of a line
 *
 * @param[in] line Line to read from
 * @param[in] number_of_entries Desired number of entries to read from,
 * starting at the beginning of the line
 *
 * @return The first desired entries of a line in a vector
 */
template <typename T>
const std::vector<T> get_multiple_entries(
    const std::string &line,
    const size_t &number_of_entries) {
  if (number_of_entries == 0)
    return {0};

  std::stringstream sstr(line);
  std::vector<T> output;

  for (size_t entry_num = 0; entry_num < number_of_entries; ++entry_num) {
    T tmp;
    sstr >> tmp;
    output.push_back(tmp);
  }

  return output;
}

/**
 * Import all entries of a line
 *
 * @param[in] line Line to read from
 *            the beginning of the line
 *
 * @return Vector of all line entries
 */
template <typename T>
const std::vector<T> get_all_entries(const std::string &line) {

  std::stringstream sstr(line);
  std::vector<T> output;

  T tmp;
  while (sstr >> tmp) {
    output.push_back(tmp);
  }

  return output;
}

/**
 * Import next line and extract the n-th entry in that line
 *
 * @param[in] file File to read from
 * @param[in] entry_number
 *            n-th entry of the next line in the given file, first line
 * entry is imported by default
 *
 * @return $n$-th entry of the next line
 */
template <typename T>
const T
get_next_line_entry(std::ifstream &file, const size_t &entry_number = 1) {
  std::string line;
  std::getline(file, line);

  return get_entry<T>(line, entry_number);
}

/**
 * Import next line and extract all entries in that line
 *
 * @param[in] file File to read from
 *
 * @return All entries of the next line in a vector
 */
template <typename T>
const std::vector<T> get_next_line_entries(std::ifstream &file) {
  std::string line;
  std::getline(file, line);

  return get_all_entries<T>(line);
}

/**
 * Import next line and extract all entries in that line, while
 * incrementing a given line number
 *
 * @param[in] file File to read from
 * @param[in] line_number Line number which will be incremented by 1
 *
 * @return All entries of the next line in a vector
 */
template <typename T>
const std::vector<T>
get_next_line_entries(std::ifstream &file, size_t &line_number) {
  std::string line;
  std::getline(file, line);
  ++line_number;

  return get_all_entries<T>(line);
}

/**
 * Skip a line in a file
 *
 * @param[in] file Input file
 * @param[in] line_skip Number of lines, which shall be skipped
 */
inline static void
skip_lines(std::ifstream &file, const size_t line_skip) {

  if (line_skip == 0)
    return;

  std::string line;
  for (size_t n = 0; n < line_skip; ++n) {
    std::getline(file, line);
  }
}

/**
 * Skip a line in a file
 *
 * @param[in] file Input file
 * @param[in] line_skip Number of lines, which shall be skipped
 * @param[in] line_number current line number, which will be incremented by
 *            the number of skipped lines;
 */
inline static void skip_lines(
    std::ifstream &file,
    const size_t line_skip,
    size_t &line_number) {

  if (line_skip == 0)
    return;

  line_number += line_skip;

  std::string line;
  for (size_t n = 0; n < line_skip; ++n) {
    std::getline(file, line);
  }
}
} // namespace Import
#endif
