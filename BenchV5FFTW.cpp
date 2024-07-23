#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <chrono>
#include <cerrno>
#include <cstring> // Include this header for strerror
#include <algorithm>
#include <sys/resource.h> // For getrusage
#include <fftw3.h>

const double PI = std::acos(-1);

using Complex = std::complex<double>;
using CArray = std::vector<Complex>;

// Function to read data from a file
void read_data(const std::string &filename, std::vector<double> &data) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening file: " << strerror(errno) << std::endl;
        return;
    }
    double value;
    while (infile >> value) {
        data.push_back(value);
    }
}

// Function to write data to a file
void write_data(const std::string &filename, const std::vector<std::complex<double>> &data) {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error opening file: " << strerror(errno) << std::endl;
        return;
    }
    for (const auto &value : data) {
        outfile << value.real() << " " << value.imag() << "\n";
    }
}

// Function to perform FFT and calculate performance metrics
void performFFT(fftw_complex *in, fftw_complex *out, int N, size_t numRuns, double corePower, bool saveOutput, bool savePerformance, std::ofstream &performanceFile) {
    double totalTime = 0;
    struct rusage usage;
    long totalMemoryUsage = 0;

    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (size_t run = 0; run < numRuns; ++run) {
        // Reset input array to original state
        for (int i = 0; i < N; ++i) {
            in[i][0] = in[i][0]; // Real part
            in[i][1] = 0.0;      // Imaginary part
        }

        auto start = std::chrono::high_resolution_clock::now();
        fftw_execute(plan);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        totalTime += elapsed.count();

        if (getrusage(RUSAGE_SELF, &usage) == 0) {
            totalMemoryUsage += usage.ru_maxrss;
        } else {
            std::cerr << "Error getting memory usage. Error code: " << errno << std::endl;
        }
    }

    double avgTime = totalTime / numRuns;
    long avgMemoryUsage = totalMemoryUsage / numRuns;

    // Estimate energy consumption
    double energyConsumptionJoules = corePower * avgTime;
    double energyConsumptionWh = energyConsumptionJoules / 3600.0;

    // Output performance metrics
    if (savePerformance) {
        performanceFile << "FFT size: " << N << "\n";
        performanceFile << "Average FFT time: " << avgTime << " seconds.\n";
        performanceFile << "Average peak memory usage: " << avgMemoryUsage << " kilobytes\n";
        performanceFile << "Estimated energy consumption per FFT run: " << energyConsumptionJoules << " joules / " << energyConsumptionWh << " watt-hours.\n\n";
    }

    // Print performance metrics to terminal
    std::cout << "FFT size: " << N << "\n";
    std::cout << "Average FFT time: " << avgTime << " seconds.\n";
    std::cout << "Average peak memory usage: " << avgMemoryUsage << " kilobytes\n";
    std::cout << "Estimated energy consumption per FFT run: " << energyConsumptionJoules << " joules / " << energyConsumptionWh << " watt-hours.\n";

    // Save FFT output to a text file only if saveOutput is true
    if (saveOutput) {
        std::ofstream outfile("fft_output.txt");
        if (!outfile) {
            std::cerr << "Error opening output file.\n";
            return;
        }

        for (int i = 0; i < N; ++i) {
            outfile << out[i][0] << " " << out[i][1] << "\n";
        }
        outfile.close();

        std::cout << "FFT output saved to fft_output.txt\n";
    }

    fftw_destroy_plan(plan);
}

int main() {
    char choice;
    std::cout << "Do you want to perform FFT on sample data (sinewave)? (y/n): ";
    std::cin >> choice;

    std::vector<double> vibrationData;
    size_t N = 0;

    if (choice == 'y' || choice == 'Y') {
        char variableChoice;
        std::cout << "Do you want to perform FFT on variable data sizes? (y/n): ";
        std::cin >> variableChoice;

        std::cout << "Enter the number of times to perform FFT: ";
        size_t numRuns;
        std::cin >> numRuns;

        std::cout << "Enter the number of processor cores: ";
        int numCores;
        std::cin >> numCores;

        std::cout << "Enter the maximum power consumption of the chip in watts: ";
        double maxChipPower;
        std::cin >> maxChipPower;

        // Calculate power consumption per core
        double corePower = maxChipPower / numCores;

        std::ofstream performanceFile("performance_measures.txt");

        if (variableChoice == 'y' || variableChoice == 'Y') {
            size_t maxPoints;
            std::cout << "Enter the maximum number of data points for variable data sizes: ";
            std::cin >> maxPoints;

            for (size_t points = 1024; points <= maxPoints; points *= 2) {
                vibrationData.resize(points);
                double frequency = 5.0; // Frequency of the sine wave
                double samplingRate = 100.0; // Sampling rate
                for (size_t i = 0; i < points; ++i) {
                    vibrationData[i] = sin(2 * PI * frequency * i / samplingRate);
                }

                N = points;
                fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
                fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

                // Fill the input array with the data
                for (int i = 0; i < N; ++i) {
                    in[i][0] = vibrationData[i]; // Real part
                    in[i][1] = 0.0;              // Imaginary part
                }

                std::cout << "Performing FFT on data of size " << N << " for " << numRuns << " times\n";
                performFFT(in, out, N, numRuns, corePower, false, true, performanceFile);

                fftw_free(in);
                fftw_free(out);
            }

            performanceFile.close();
            std::cout << "Performance measures saved to performance_measures.txt\n";
        } else {
            std::cout << "Enter the number of data points for fixed data size: ";
            std::cin >> N;

            vibrationData.resize(N);
            double frequency = 5.0; // Frequency of the sine wave
            double samplingRate = 100.0; // Sampling rate
            for (size_t i = 0; i < N; ++i) {
                vibrationData[i] = sin(2 * PI * frequency * i / samplingRate);
            }

            fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
            fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

            // Fill the input array with the data
            for (int i = 0; i < N; ++i) {
                in[i][0] = vibrationData[i]; // Real part
                in[i][1] = 0.0;              // Imaginary part
            }

            std::cout << "Performing FFT on data of size " << N << " for " << numRuns << " times\n";
            performFFT(in, out, N, numRuns, corePower, true, false, performanceFile);

            fftw_free(in);
            fftw_free(out);
        }
    } else {
        std::ifstream infile("vibration_data.txt");
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

        fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

        // Fill the input array with the data
        for (int i = 0; i < N; ++i) {
            in[i][0] = vibrationData[i]; // Real part
            in[i][1] = 0.0;              // Imaginary part
        }

        std::cout << "Enter the number of times to perform FFT: ";
        size_t numRuns;
        std::cin >> numRuns;

        std::cout << "Enter the number of processor cores: ";
        int numCores;
        std::cin >> numCores;

        std::cout << "Enter the maximum power consumption of the chip in watts: ";
        double maxChipPower;
        std::cin >> maxChipPower;

        // Calculate power consumption per core
        double corePower = maxChipPower / numCores;

        std::cout << "Performing FFT on vibration data of size " << N << " for " << numRuns << " times\n";
        std::ofstream performanceFile("performance_measures.txt");
        performFFT(in, out, N, numRuns, corePower, true, true, performanceFile);
        performanceFile.close();

        fftw_free(in);
        fftw_free(out);
    }

    return 0;
}
