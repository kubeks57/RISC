#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <cerrno>

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
    // Define sampling rate
    double Fs = 1000.0; // Example: 1000 Hz, replace with your actual sampling rate

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

    // Perform FFT on the data
    fft(data);

    // Save the FFT output to a text file
    std::ofstream outfile("fft_PC_output.txt");
    if (!outfile.is_open()) {
        std::cerr << "Error opening output file. Error code: " << errno << std::endl;
        return 1;
    }

    for (size_t i = 0; i < nextPowerOf2; ++i) {
        outfile << data[i].real() << "\t" << data[i].imag() << "\n";
    }
    outfile.close();

    std::cout << "FFT output saved to 'fft_PC_output.txt'.\n";

    return 0;
}
