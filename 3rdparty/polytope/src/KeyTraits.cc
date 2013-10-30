//-----------------------------------------------------------------------------//
// KeyTraits
//
// Encapsulate how we think about keys for the space filling curves.
// This is specialized assuming we're working with 64 bit unsigned long long.
//----------------------------------------------------------------------------//
#include "KeyTraits.hh"

namespace polytope {

const uint32_t KeyTraits::numbits = 64U;
const uint32_t KeyTraits::numbits1d = 20U;
const KeyTraits::Key KeyTraits::zero = 0ULL;
const KeyTraits::Key KeyTraits::one = 1ULL;
const KeyTraits::Key KeyTraits::two = 2ULL;
const KeyTraits::Key KeyTraits::maxKey1d = 2097152ULL;
const KeyTraits::Key KeyTraits::maxKey = 9223372036854775808ULL;

}
