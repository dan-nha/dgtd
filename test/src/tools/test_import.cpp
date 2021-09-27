#include "../../../src/tools/import.h"

#include <boost/test/unit_test.hpp>
#include <fstream>
#include <string>
#include <vector>

BOOST_AUTO_TEST_SUITE(tools);

BOOST_AUTO_TEST_CASE(get_line_entries) {
  using namespace Import;

  std::string eminem("Hi! My name is Slim Shady.");

  BOOST_TEST(get_entry<std::string>(eminem) == "Hi!");
  BOOST_TEST(get_entry<std::string>(eminem, 5) == "Slim");

  std::vector<std::string> partial_line(
      get_multiple_entries<std::string>(eminem, 4));
  BOOST_TEST(partial_line[0] == "Hi!");
  BOOST_TEST(partial_line[1] == "My");
  BOOST_TEST(partial_line[2] == "name");
  BOOST_TEST(partial_line[3] == "is");

  std::vector<std::string> total_line(
      get_all_entries<std::string>(eminem));
  BOOST_TEST(total_line[0] == "Hi!");
  BOOST_TEST(total_line[1] == "My");
  BOOST_TEST(total_line[2] == "name");
  BOOST_TEST(total_line[3] == "is");
  BOOST_TEST(total_line[4] == "Slim");
  BOOST_TEST(total_line[5] == "Shady.");
}

BOOST_AUTO_TEST_CASE(get_line_entries_from_file) {
  using namespace Import;

  const std::string root_dir(DGTD_ROOT);
  std::string filename(root_dir + "/test/src/tools/test.txt");
  std::ifstream file(filename.c_str());

  BOOST_TEST(get_next_line_entry<std::string>(file) == "Hi!");

  std::vector<std::string> lyrics(
      get_next_line_entries<std::string>(file));
  BOOST_TEST(lyrics[2] == "is");

  skip_lines(file, 2);
  BOOST_TEST(get_next_line_entries<std::string>(file).back() == "(huh?)");

  size_t line_number(5);
  skip_lines(file, 3, line_number);
  BOOST_TEST(get_next_line_entries<std::string>(file).back() == "(what?)");
  BOOST_TEST(line_number == 8);
}

BOOST_AUTO_TEST_SUITE_END();
