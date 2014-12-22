#include "mathutils.h"
namespace ham {
// ----------------------------------------------------------------------------------------
// add two numbers, treating -INFINITY as zero, i.e. calculates log a*b = log a + log b, i.e. a *and* b
double AddWithMinusInfinities(double first, double second) {
  if(first == -INFINITY || second == -INFINITY)
    return -INFINITY;
  else
    return first + second;
}

}
