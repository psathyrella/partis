// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"


namespace igsw {

TEST_CASE("Trivial pass", "[trivial]") {
  REQUIRE(1 == 1);
}

TEST_CASE("Trivial fail", "[trivial]") {
  REQUIRE(1 == 2);
}

}
