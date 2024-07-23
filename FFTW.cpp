#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <fftw3.h>

void read_data(const std::string &filename, std::vector<double> &data) {
    std::ifstream infile(filename);
    double value;
    while (infile >> value) {
        data.push_back(value);
    }
}

void write_data(const std::string &filename, const std::vector<std::complex<double>> &data) {
    std::ofstream outfile(filename);
    for (const auto &value : data) {
        outfile << value.real() << " " << value.imag() << "\n";
    }
}

int main() {
    std::string input_filename = "vibration_data.txt";
    std::string output_filename = "FFTW_Out.txt";

    // Read input data
    std::vector<double> input_data;
    read_data(input_filename, input_data);

    int N = input_data.size();

    // Allocate input and output arrays for FFTW
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Fill the input array with the data
    for (int i = 0; i < N; ++i) {
        in[i][0] = input_data[i]; // Real part
        in[i][1] = 0.0;           // Imaginary part
    }

    // Plan and execute the FFT
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Store the results in a vector of complex numbers
    std::vector<std::complex<double>> output_data(N);
    for (int i = 0; i < N; ++i) {
        output_data[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    // Write the output data to a file
    write_data(output_filename, output_data);

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    std::cout << "FFT computation complete. Results saved to " << output_filename << std::endl;
    return 0;
}
