#include <array>
#include <cstdint>
#include <vector>
#include <x86intrin.h>

template <std::size_t N, std::size_t... J>
std::size_t get_index_morton_helper(const std::array<std::size_t, N> &i,
                                    const std::array<std::size_t, N> &m,
                                    std::index_sequence<J...>) {
  return (_pdep_u64(i[J], m[J]) | ...);
}

template <std::size_t N>
std::size_t get_index_morton(const std::vector<std::array<std::size_t, N>> i,
                             const std::array<std::size_t, N> m) {
  for (std::size_t j = 0; j < 2048; ++j) {
    __asm volatile("# LLVM-MCA-BEGIN Morton%0" ::[N] "i"(N));

    std::size_t r =
        get_index_morton_helper(i[j], m, std::make_index_sequence<N>());

    __asm volatile("# LLVM-MCA-END Morton%0" ::[N] "i"(N), "r"(r));
  }

  return 0;
}

template <std::size_t N, std::size_t I>
std::size_t get_index_canonical_helper(const std::array<std::size_t, N> &i,
                                       const std::array<std::size_t, N> &s) {
  if constexpr (I == N - 1) {
    return i[I];
  } else {
    return i[I] + s[I] * get_index_canonical_helper<N, I + 1>(i, s);
  }
}

template <std::size_t N>
std::size_t get_index_canonical(const std::vector<std::array<std::size_t, N>> i,
                                const std::array<std::size_t, N> s) {
  for (std::size_t j = 0; j < 2048; ++j) {
    __asm volatile("# LLVM-MCA-BEGIN Canonical%0" ::[N] "i"(N));

    std::size_t r = get_index_canonical_helper<N, 0>(i[j], s);

    __asm volatile("# LLVM-MCA-END Canonical%0" ::[N] "i"(N), "r"(r));
  }

  return 0;
}

#define MORTON(n)                                                              \
  template std::size_t get_index_morton<n>(                                    \
      const std::vector<std::array<std::size_t, n>>,                           \
      const std::array<std::size_t, n>)
#define CANONICAL(n)                                                           \
  template std::size_t get_index_canonical<n>(                                 \
      const std::vector<std::array<std::size_t, n>>,                           \
      const std::array<std::size_t, n>)

MORTON(1);
MORTON(2);
MORTON(3);
MORTON(4);
MORTON(5);
MORTON(6);
MORTON(7);
MORTON(8);
MORTON(9);
MORTON(10);

CANONICAL(1);
CANONICAL(2);
CANONICAL(3);
CANONICAL(4);
CANONICAL(5);
CANONICAL(6);
CANONICAL(7);
CANONICAL(8);
CANONICAL(9);
CANONICAL(10);
