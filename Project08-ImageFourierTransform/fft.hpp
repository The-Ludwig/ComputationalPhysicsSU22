#pragma once

#ifndef NDEBUG
#include <iostream>
#endif

#include <complex>
#include <vector>

using std::sqrt, std::exp, std::abs, std::complex, std::size_t;
using std::vector;
constexpr std::complex<double> i(0, 1);
constexpr double pi = 3.14159265358979323846;
constexpr double isqrt2 = 0.7071067811865475;

size_t bitReverse(size_t x, int log2n) {
  int n = 0;
  for (int i = 0; i < log2n; i++) {
    n <<= 1;
    n |= (x & 1);
    x >>= 1;
  }
  return n;
}

/**
 * @brief Fourier transform for arbitrary input length
 *
 */
template <class SrcIter_t, class DstIter_t>
void dft(const SrcIter_t src, DstIter_t dest, size_t len) {
  typedef typename std::iterator_traits<DstIter_t>::value_type cmp_dest;

  cmp_dest i = cmp_dest(0, 1);

  for (size_t k = 0; k < len; k++) {
    cmp_dest bf = exp(-i * 2. * pi * double(k) / double(len));
    cmp_dest f = cmp_dest(1, 0);
    dest[k] = 0;
    for (size_t n = 0; n < len; n++) {
      dest[k] += cmp_dest(src[n]) * f;
      f *= bf;
    }
    dest[k] /= sqrt(double(len));
  }
}

template <class SrcIter_t, class DstIter_t>
void fft_power_of_two(const SrcIter_t src, DstIter_t dest, size_t len) {
  typedef typename std::iterator_traits<DstIter_t>::value_type cmp_dest;

  size_t log2n = 0;
  while ((len >> ++log2n) != 0) {
  };
  log2n--;

  // reorder (bit_reverse)
  size_t m = 0;
  for (size_t k = 0; k < len; k++)
    dest[k] = cmp_dest(src[bitReverse(m, log2n)]);

  for (size_t s = 1; s <= log2n; s++) {
    size_t m = 1 << s;
    size_t m2 = m >> 1;
    cmp_dest w(1, 0);
    cmp_dest wm = exp(-i * (pi / m2));
    for (size_t j = 0; j < m2; ++j) {
      for (size_t k = j; k < len; k += m) {
        cmp_dest a = dest[k];
        cmp_dest b = w * dest[k + m2];
        dest[k] = (a + b) * isqrt2;
        dest[k + m2] = (a - b) * isqrt2;
      }
      w *= wm;
    }
  }
}

size_t new_index(size_t idx, size_t len, size_t npf2) {
  size_t ni = idx;

  for (size_t k = 0; k < npf2; k++) {
    size_t base = (ni / (len >> k)) * (len >> k);
    ni = (len >> (k + 1)) * (ni % 2) + base + (ni - base) / 2;
  }

  return ni;
}

template <class SrcIter_t, class DstIter_t>
void fft(const SrcIter_t src, DstIter_t dest, size_t len) {
  typedef typename std::iterator_traits<DstIter_t>::value_type cmp_dest;

  size_t log2n = 0;
  while (((len >> ++log2n) & 1) == 0) {
  };
  log2n--;

  // reorder
  for (size_t k = 0; k < len; k++) {
    dest[new_index(k, len, log2n)] = cmp_dest(src[k]);
  }

  // fourier traffos of inner parts, which are not divisible by two
  size_t inner_len = len / (1 << log2n);
  // buffer
  cmp_dest* buf = new cmp_dest[inner_len];
  for (size_t p = 0; p < (1u << log2n); p++) {
    dft(&(dest[inner_len * p]), buf, inner_len);
    for (size_t k = 0; k < inner_len; k++) dest[inner_len * p + k] = buf[k];
  }
  delete[] buf;

  for (size_t s = 1; s <= log2n; s++) {
    size_t m = inner_len * (1u << s);
    size_t m2 = m >> 1;
    cmp_dest w(1, 0);
    cmp_dest wm = exp(-i * (pi / m2));
    for (size_t j = 0; j < m2; ++j) {
      for (size_t k = j; k < len; k += m) {
        cmp_dest a = dest[k];
        cmp_dest b = w * dest[k + m2];
        dest[k] = (a + b) * isqrt2;
        dest[k + m2] = (a - b) * isqrt2;
      }
      w *= wm;
    }
  }
}

template <class SrcIter_t, class DstIter_t>
void ifft(const SrcIter_t src, DstIter_t dest, size_t len) {
  typedef typename std::iterator_traits<DstIter_t>::value_type cmp_dest;
  cmp_dest* buf = new cmp_dest[len];
  fft(src, buf, len);
  dest[0] = buf[0];
  for (size_t k = 1; k < len; k++) dest[k] = buf[len - k];
  delete[] buf;
}

template <class SrcIter_t, class DstIter_t>
void idft(const SrcIter_t src, DstIter_t dest, size_t len) {
  typedef typename std::iterator_traits<DstIter_t>::value_type cmp_dest;
  cmp_dest* buf = new cmp_dest[len];
  dft(src, buf, len);
  dest[0] = buf[0];
  for (size_t k = 1; k < len; k++) dest[k] = buf[len - k];
  delete[] buf;
}

/**
 * @brief Discrete cosine transform for arbitrary input length
 *
 */
template <class SrcIter_t, class DstIter_t>
void dct(const SrcIter_t src, DstIter_t dest, size_t len) {
  typedef typename std::iterator_traits<DstIter_t>::value_type cmp_dest;

  // reorder
  for (size_t k = 0; k < len / 2 - 1; k++) {
    dest[k] = src[2 * k];
    dest[len - 1 - k] = src[2 * k + 1];
  }
  if (len % 2 == 1) dest[len / 2] = src[len - 1];

  std::complex<cmp_dest>* buf = new std::complex<cmp_dest>[len];
  ifft(dest, buf, len);

  dest[0] = 2 * isqrt2 * std::real(buf[0]);

  std::complex<cmp_dest> w(1, 0);
  std::complex<cmp_dest> wm = exp(i * pi / (2. * len));
  for (size_t k = 1; k < len; k++) {
    w *= wm;
    dest[k] = 2 * std::real(w * buf[k]);
  }
  delete[] buf;
}

#ifndef NDEBUG
namespace FFT {

void simple() {
  std::cout << "simple power of two test" << std::endl;
  double a[] = {0, 1, 2, 3, 2, 1, 0, 0};
  complex<double> b[8];
  fft_power_of_two(a, b, 8);
  for (int i = 0; i < 8; ++i) std::cout << b[i] << "\n";
  std::cout << std::endl;
  dft(a, b, 8);
  for (int i = 0; i < 8; ++i) std::cout << b[i] << "\n";
  std::cout << std::endl;
  fft(a, b, 8);
  for (int i = 0; i < 8; ++i) std::cout << b[i] << "\n";
}

void prime() {
  std::cout << "prime test" << std::endl;
  double a[] = {0, 1, 2, 3, 2, 1, 0};
  complex<double> b[7];

  std::cout << std::endl;
  dft(a, b, 7);
  for (int i = 0; i < 7; ++i) std::cout << b[i] << "\n";
  std::cout << std::endl;
  fft(a, b, 7);
  for (int i = 0; i < 7; ++i) std::cout << b[i] << "\n";
}

void interesting() {
  std::cout << "interesting test" << std::endl;
  double a[] = {0, 1, 2, 3, 2, 1, 0, 5, 9, 1, 1, -3};
  complex<double> b[12];

  std::cout << std::endl;
  dft(a, b, 12);
  for (int i = 0; i < 12; ++i) std::cout << b[i] << "\n";
  std::cout << std::endl;
  fft(a, b, 12);
  for (int i = 0; i < 12; ++i) std::cout << b[i] << "\n";
}

void inverse() {
  // Testing the inverse FFT
  std::cout << "inverse test" << std::endl;
  typedef complex<double> cx;
  cx a[] = {cx(0), cx(1), cx(2), cx(3), cx(2), cx(1),
            cx(0), cx(5), cx(9), cx(1), cx(1), cx(-3)};
  cx b[12];

  std::cout << std::endl;
  fft(a, b, 12);
  for (int i = 0; i < 12; ++i) std::cout << b[i] << "\n";

  std::cout << std::endl;
  for (int i = 0; i < 12; ++i) std::cout << a[i] << "\n";
  std::cout << std::endl;
  ifft(b, a, 12);
  for (int i = 0; i < 12; ++i) std::cout << a[i] << "\n";
}

void tests() {
  simple();
  prime();
  interesting();
  inverse();
}
};  // namespace FFT
#endif