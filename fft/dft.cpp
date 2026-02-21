#include <vector>
#include <complex>
#include <cmath>
#include <numbers>
#include <iostream>
#include <string>

using cd = std::complex<double>;
constexpr double PI = std::numbers::pi;

std::vector<cd> dft(const std::vector<cd>& x) {
	int n = static_cast<int>(x.size());
	std::vector<cd> ft(n, cd(0.0, 0.0));
	for (int k = 0; k < n; k++) {
		cd sum{0.0, 0.0};
		for (int m = 0; m < n; m++) {
			const cd angle = cd(0.0, -2.0 * PI * k * m / n);
			sum += x[m] * std::exp(angle);
		}
		ft[k] = sum;
	}
	return ft;
}

void fft_rec(cd* input, cd* output, int N, int stride) {
    if (N == 1) {
        output[0] = input[0];
    }
    else {
        fft_rec(input, output, N / 2, 2 * stride);
        fft_rec(input + stride, output + N / 2, N / 2, 2 * stride);
        for (int k = 0; k < (N / 2); k++) {
            // need to implement.
            cd p = output[k];
            cd angle = cd(0.0, -2.0 * PI * k / N);
            cd q = std::exp(angle) * output[k + N / 2];
            output[k] = p + q;
            output[k + N / 2] = p - q;
        }
    }
}

// n must be a power of 2
std::vector<cd> ctfft(const std::vector<cd>& input){
    std::vector<cd> input_ = input;
    std::vector<cd> output(input.size());
    fft_rec(input_.data(), output.data(), input_.size(), 1);
    return output;
}

void print_vec(const std::string& label, const std::vector<cd>& v) {
    std::cout << label << "\n";
    for (const auto& x : v) {
        std::cout << x << "\n";
    }
    std::cout << "\n";
}

void test_case(const std::string& name,
    const std::vector<cd>& input,
    const std::string& expected) {
    std::cout << name << "\n";
    auto a = dft(input);
    auto b = ctfft(input);
    print_vec("dft:", a);
    print_vec("ctfft:", b);
    std::cout << "Expected: " << expected << "\n\n";
}

int main() {
    test_case(
        "Impulse test:",
        { cd(1,0), cd(0,0), cd(0,0), cd(0,0) },
        "all (1,0)"
    );

    test_case(
        "Constant signal test:",
        { cd(1,0), cd(1,0), cd(1,0), cd(1,0) },
        "(4,0) then three (0,0)"
    );

    test_case(
        "Single-frequency test:",
        { cd(1,0), cd(0,0), cd(-1,0), cd(0,0) },
        "(0,0), (2,0), (0,0), (2,0)"
    );

    return 0;
}