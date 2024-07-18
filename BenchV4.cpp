#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <sys/resource.h> // For getrusage

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

// Function to perform FFT for variable data sizes
void performVariableSizeFFT(const std::vector<double>& vibrationData, size_t maxDataPoints, size_t numRuns, int numCores, double maxChipPower) {
    std::ofstream performanceFile("fft_performance_results.txt");
    if (!performanceFile) {
        std::cerr << "Error opening performance file.\n";
        return;
    }

    for (size_t dataSize = 1024; dataSize <= maxDataPoints; dataSize *= 2) {
        std::cout << "Performing FFT on data of size " << dataSize << " for " << numRuns << " times\n";

        size_t nextPowerOf2 = std::pow(2, std::ceil(std::log2(dataSize)));
        std::vector<Complex> expTable(nextPowerOf2 / 2);
        for (size_t i = 0; i < nextPowerOf2 / 2; ++i) {
            expTable[i] = std::exp(Complex(0, -2 * PI * i / nextPowerOf2));
        }

        std::vector<unsigned int> bitReversal(nextPowerOf2);
        precomputeBitReversal(bitReversal, std::log2(nextPowerOf2));

        double totalTime = 0;
        size_t totalComplexAdds = 0;
        size_t totalComplexMults = 0;

        for (size_t run = 0; run < numRuns; ++run) {
            CArray data(nextPowerOf2);
            for (size_t i = 0; i < nextPowerOf2; ++i) {
                data[i] = Complex(vibrationData[i % dataSize], 0);
            }

            auto start = std::chrono::high_resolution_clock::now();
            fft(data, expTable, bitReversal);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            totalTime += elapsed.count();

            totalComplexAdds += complex_adds;
            totalComplexMults += complex_mults;
            complex_adds = 0; // Reset for the next run
            complex_mults = 0; // Reset for the next run
        }

        double avgTime = totalTime / numRuns;
        size_t avgComplexAdds = totalComplexAdds / numRuns;
        size_t avgComplexMults = totalComplexMults / numRuns;

        performanceFile << "Data Size: " << dataSize << "\n";
        performanceFile << "Average FFT time: " << avgTime << " seconds.\n";
        performanceFile << "Estimated number of complex multiplications on data: " << avgComplexMults << '\n';
        performanceFile << "Estimated number of complex additions on data: " << avgComplexAdds << '\n';
        performanceFile << "\n";
    }

    performanceFile.close();
    std::cout << "Performance results saved to fft_performance_results.txt\n";
}

int main() {
    char choice;
    std::cout << "Do you want to perform FFT on variable number of data points? (y/n): ";
    std::cin >> choice;

    std::vector<double> vibrationData;
    size_t N;

    if (choice == 'y' || choice == 'Y') {
        size_t maxDataPoints;
        std::cout << "Enter the maximum number of data points for variable FFT: ";
        std::cin >> maxDataPoints;
        
        vibrationData.resize(maxDataPoints);
        double frequency = 5.0; // Frequency of the sine wave
        double samplingRate = 100.0; // Sampling rate
        for (size_t i = 0; i < maxDataPoints; ++i) {
            vibrationData[i] = sin(2 * PI * frequency * i / samplingRate);
        }

        size_t numRuns;
        std::cout << "Enter the number of times to perform FFT: ";
        std::cin >> numRuns;
        int numCores;
        std::cout << "Enter the number of processor cores: ";
        std::cin >> numCores;
        double maxChipPower;
        std::cout << "Enter the maximum power consumption of the chip in watts: ";
        std::cin >> maxChipPower;
        
        performVariableSizeFFT(vibrationData, maxDataPoints, numRuns, numCores, maxChipPower);
        return 0; // Exit after performing variable size FFT
    } else {
        std::cout << "Do you want to perform FFT on data from a file? (y/n): ";
        std::cin >> choice;

        if (choice == 'y' || choice == 'Y') {
            std::string filename;
            std::cout << "Enter the filename: ";
            std::cin >> filename;

            std::ifstream infile(filename);
            if (!infile) {
                std::cerr << "Error opening input file.\n";
                return 1;
            }

            double value;
            while (infile >> value) {
                vibrationData.push_back(value);
            }
            infile.close();
            
            N = vibrationData.size();
        } else {
            std::cout << "Enter the number of data points: ";
            std::cin >> N;
            
            vibrationData.resize(N);
            double frequency = 5.0; // Frequency of the sine wave
            double samplingRate = 100.0; // Sampling rate
            for (size_t i = 0; i < N; ++i) {
                vibrationData[i] = sin(2 * PI * frequency * i / samplingRate);
            }
        }
    }

    // Proceed with fixed size FFT as originally implemented
    size_t nextPowerOf2 = std::pow(2, std::ceil(std::log2(N)));
    vibrationData.resize(nextPowerOf2, 0.0); // Pad with zeros

    // Convert vibration data to complex array
    CArray data(nextPowerOf2);
    for (size_t i = 0; i < nextPowerOf2; ++i) {
        data[i] = Complex(vibrationData[i], 0);
    }

    size_t numRuns;
    std::cout << "Enter the number of times to perform FFT: ";
    std::cin >> numRuns;

    int numCores;
    std::cout << "Enter the number of processor cores: ";
    std::cin >> numCores;

    double maxChipPower;
    std::cout << "Enter the maximum power consumption of the chip in watts: ";
    std::cin >> maxChipPower;

    // Calculate power consumption per core
    double corePower = maxChipPower / numCores;

    // Precompute complex exponentials and bit-reversal indices
    std::vector<Complex> expTable(nextPowerOf2 / 2);
    for (size_t i = 0; i < nextPowerOf2 / 2; ++i) {
        expTable[i] = std::exp(Complex(0, -2 * PI * i / nextPowerOf2));
    }

    std::vector<unsigned int> bitReversal(nextPowerOf2);
    precomputeBitReversal(bitReversal, std::log2(nextPowerOf2));

    // Measure time for FFT and calculate average time
    double totalTime = 0;
    size_t totalComplexAdds = 0;
    size_t totalComplexMults = 0;
    struct rusage usage;
    long totalMemoryUsage = 0;

    for (size_t run = 0; run < numRuns; ++run) {
        // Reset data to original state
        for (size_t i = 0; i < nextPowerOf2; ++i) {
            data[i] = Complex(vibrationData[i], 0);
        }

        auto start = std::chrono::high_resolution_clock::now();
        fft(data, expTable, bitReversal);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        totalTime += elapsed.count();

        totalComplexAdds += complex_adds;
        totalComplexMults += complex_mults;
        complex_adds = 0; // Reset for the next run
        complex_mults = 0; // Reset for the next run

        if (getrusage(RUSAGE_SELF, &usage) == 0) {
            totalMemoryUsage += usage.ru_maxrss;
        } else {
            std::cerr << "Error getting memory usage. Error code: " << errno << std::endl;
        }
    }

    double avgTime = totalTime / numRuns;
    size_t avgComplexAdds = totalComplexAdds / numRuns;
    size_t avgComplexMults = totalComplexMults / numRuns;
    long avgMemoryUsage = totalMemoryUsage / numRuns;

    std::cout << "Average FFT time: " << avgTime << " seconds.\n";
    std::cout << "Average peak memory usage: " << avgMemoryUsage << " kilobytes\n";
    std::cout << "Estimated number of complex multiplications on data: " << avgComplexMults << '\n';
    std::cout << "Estimated number of complex additions on data: " << avgComplexAdds << '\n';

    // Calculate total real operations
    size_t totalRealMults = avgComplexMults * 4;
    size_t totalRealAdds = avgComplexMults * 2 + avgComplexAdds * 2;
    std::cout << "Estimated total number of real operations on data: " << (totalRealMults + totalRealAdds) << '\n';

    // Estimate energy consumption
    double energyConsumptionJoules = corePower * avgTime;
    double energyConsumptionWh = energyConsumptionJoules / 3600.0;
    
    std::cout << "Estimated energy consumption per FFT run: " << energyConsumptionJoules << " joules/ " << energyConsumptionWh << " watt-hours.\n";

    // Save FFT output to a text file
    std::ofstream outfile("fft_output.txt");
    if (!outfile) {
        std::cerr << "Error opening output file.\n";
        return 1;
    }

    for (const auto& val : data) {
        outfile << val.real() << " " << val.imag() << "\n";
    }
    outfile.close();

    std::cout << "FFT output saved to fft_output.txt\n";

    return 0;
}
