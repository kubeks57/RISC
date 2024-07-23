#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

// Function to read data from a text file
std::vector<double> readData(const std::string& filename) {
    std::vector<double> data;
    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Error opening input file." << std::endl;
        return data;
    }

    double value;
    while (inputFile >> value) {
        data.push_back(value);
    }

    inputFile.close();
    return data;
}

// Function to write data to a text file
void writeData(const std::string& filename, const std::vector<std::complex<double>>& data) {
    std::ofstream outputFile(filename);
    if (!outputFile) {
        std::cerr << "Error opening output file." << std::endl;
        return;
    }

    for (const auto& value : data) {
        outputFile << value.real() << " " << value.imag() << std::endl;
    }

    outputFile.close();
}

// Function to compute the DFT
std::vector<std::complex<double>> computeDFT(const std::vector<double>& data) {
    int N = data.size();
    std::vector<std::complex<double>> dft(N);

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double angle = 2 * M_PI * k * n / N;
            std::complex<double> expTerm(std::cos(angle), -std::sin(angle));
            sum += data[n] * expTerm;
        }
        dft[k] = sum;
    }

    return dft;
}

int main() {
    std::string inputFilename = "vibration_data.txt";
    std::string outputFilename = "DFT.txt";

    // Read data from input file
    std::vector<double> data = readData(inputFilename);
    if (data.empty()) {
        return 1;
    }

    // Compute the DFT
    std::vector<std::complex<double>> dft = computeDFT(data);

    // Write DFT results to output file
    writeData(outputFilename, dft);

    std::cout << "DFT computation and output completed successfully." << std::endl;
    return 0;
}
