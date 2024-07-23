#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>

// Function to read FFT data from a file
void read_fft_data(const std::string &filename, std::vector<std::complex<double>> &data) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    double real, imag;
    while (infile >> real >> imag) {
        data.emplace_back(real, imag);
    }
}

// Function to normalize FFT data
void normalize_data(std::vector<std::complex<double>> &data) {
    double max_val = 0.0;
    for (const auto &value : data) {
        max_val = std::max(max_val, std::abs(value));
    }

    if (max_val > 0) {
        for (auto &value : data) {
            value /= max_val;
        }
    }
}

// Function to compute error metrics
void compute_error_metrics(const std::vector<std::complex<double>> &data1, const std::vector<std::complex<double>> &data2) {
    size_t N1 = data1.size();
    size_t N2 = data2.size();

    // Determine the size to compare (minimum of the two sizes)
    size_t minSize = std::min(N1, N2);

    // Calculate error metrics
    double mae = 0.0;
    double mse = 0.0;
    double max_abs_error = 0.0;

    for (size_t i = 0; i < minSize; ++i) {
        double abs_error = std::abs(data1[i] - data2[i]);
        mae += abs_error;
        mse += abs_error * abs_error;
        max_abs_error = std::max(max_abs_error, abs_error);
    }

    mae /= minSize;
    mse /= minSize;
    double rmse = std::sqrt(mse);

    std::cout << "Mean Absolute Error (MAE): " << mae << std::endl;
    std::cout << "Root Mean Squared Error (RMSE): " << rmse << std::endl;
    std::cout << "Maximum Absolute Error: " << max_abs_error << std::endl;
}

int main() {
    std::string file1, file2;
    std::cout << "Enter the filename of the first FFT result: ";
    std::cin >> file1;
    std::cout << "Enter the filename of the second FFT result: ";
    std::cin >> file2;

    std::vector<std::complex<double>> data1, data2;

    // Read FFT data from files
    read_fft_data(file1, data1);
    read_fft_data(file2, data2);

    // Normalize the data
    normalize_data(data1);
    normalize_data(data2);

    // Compute error metrics
    compute_error_metrics(data1, data2);

    return 0;
}
