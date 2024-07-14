#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <cerrno>
#include <algorithm>

const double PI = std::acos(-1);

using Complex = std::complex<double>;
using CArray = std::vector<Complex>;

// Global counters for operations
size_t complex_adds = 0;
size_t complex_mults = 0;

// Precompute bit-reversal indices for bit-reversal permutation
void precomputeBitReversal(std::vector<unsigned int>& bitReversal, int log2n) {
    const size_t N = bitReversal.size();
    for (size_t i = 0; i < N; ++i) {
        unsigned int x = i;
        unsigned int n = 0;
        for (int j = 0; j < log2n; ++j) {
            n <<= 1;
            n |= (x & 1);
            x >>= 1;
        }
        bitReversal[i] = n;
    }
}

// Iterative FFT using the Cooley-Tukey algorithm
void fft(CArray& x, const std::vector<Complex>& expTable, const std::vector<unsigned int>& bitReversal) {
    const size_t N = x.size();
    const int log2n = std::log2(N);

    // Bit-reversal permutation using precomputed table
    for (size_t i = 0; i < N; ++i) {
        size_t rev = bitReversal[i];
        if (i < rev) {
            std::swap(x[i], x[rev]);
        }
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
                
                // Increment the operation counters
                complex_mults++; // one complex multiplication
                complex_adds += 2; // two complex additions
            }
        }
    }
}

int main() {
    double Fs = 20000; // Sampling rate

    // Read vibration data from file
    std::ifstream infile("vibration_data.txt");
    if (!infile) {
        std::cerr << "Error opening input file.\n";
        return 1;
    }

    std::vector<double> vibrationData;
    double value;
    while (infile >> value) {
        vibrationData.push_back(value);
    }
    infile.close();

    // Pad data to the nearest power of 2
    size_t N = vibrationData.size();
    size_t nextPowerOf2 = std::pow(2, std::ceil(std::log2(N)));
    vibrationData.resize(nextPowerOf2, 0.0); // Pad with zeros

    // Convert vibration data to complex array
    CArray data(nextPowerOf2);
    for (size_t i = 0; i < nextPowerOf2; ++i) {
        data[i] = Complex(vibrationData[i], 0);
    }

    std::cout << "Performing FFT on data of size " << nextPowerOf2 << '\n';

    // Precompute complex exponentials and bit-reversal indices
    std::vector<Complex> expTable(nextPowerOf2 / 2);
    for (size_t i = 0; i < nextPowerOf2 / 2; ++i) {
        expTable[i] = std::exp(Complex(0, -2 * PI * i / nextPowerOf2));
    }

    std::vector<unsigned int> bitReversal(nextPowerOf2);
    precomputeBitReversal(bitReversal, std::log2(nextPowerOf2));

    // Measure time for FFT
    auto start = std::chrono::high_resolution_clock::now();
    fft(data, expTable, bitReversal);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "FFT took " << elapsed.count() << " seconds.\n";

    // Save the FFT output to a text file
    std::ofstream outfile("fft_PC_output_time.txt");
    if (!outfile.is_open()) {
        std::cerr << "Error opening output file. Error code: " << errno << std::endl;
        return 1;
    }

    for (size_t i = 0; i < nextPowerOf2; ++i) {
        outfile << data[i].real() << "\t" << data[i].imag() << "\n";
    }
    outfile.close();

    std::cout << "FFT output saved to 'fft_PC_output_time.txt'.\n";

    // Print the number of operations
    std::cout << "Number of complex multiplications: " << complex_mults << '\n';
    std::cout << "Number of complex additions: " << complex_adds << '\n';
    std::cout << "Total number of operations: " << (complex_mults * 6 + complex_adds * 2) << '\n';

    return 0;
}
