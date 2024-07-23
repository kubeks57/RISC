#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <chrono>

const double PI = std::acos(-1);

using Complex = std::complex<double>;
using CArray = std::vector<Complex>;

// Function to reverse bits for bit-reversal permutation
unsigned int reverseBits(unsigned int x, int log2n) {
    unsigned int n = 0;
    for (int i = 0; i < log2n; ++i) {
        n <<= 1;
        n |= (x & 1);
        x >>= 1;
    }
    return n;
}

// Iterative FFT using the Cooley-Tukey algorithm
void fft(CArray& x) {
    const size_t N = x.size();
    const int log2n = std::log2(N);

    // Bit-reversal permutation
    for (size_t i = 0; i < N; ++i) {
        size_t rev = reverseBits(i, log2n);
        if (i < rev) {
            std::swap(x[i], x[rev]);
        }
    }

    // Precompute complex exponentials
    std::vector<Complex> expTable(N / 2);
    for (size_t i = 0; i < N / 2; ++i) {
        expTable[i] = std::exp(Complex(0, -2 * PI * i / N));
    }

    // FFT
    for (size_t s = 1; s <= log2n; ++s) {
        size_t m = 1 << s;
        size_t m2 = m >> 1;
        for (size_t j = 0; j < m2; ++j) {
            Complex w = expTable[j * N / m];
            for (size_t k = j; k < N; k += m) {
                Complex t = w * x[k + m2];
                Complex u = x[k];
                x[k] = u + t;
                x[k + m2] = u - t;
            }
        }
    }
}

int main() {
    const size_t N = 1024; // Input size
    CArray data;
    data.reserve(N);

    // Sample data - sine wave
    for (size_t i = 0; i < N; ++i) {
        data.emplace_back(std::sin(2 * PI * i / N), 0);
    }

    std::cout << "Performing FFT on data of size " << N << '\n';

    // Benchmarking 1000 runs
    const int num_runs = 1000;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_runs; ++i) {
        CArray data_copy = data; // Create a copy for each run
        fft(data_copy);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Average elapsed time: " << (elapsed.count() / num_runs) << " seconds\n";

    return 0;
}
