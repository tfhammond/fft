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



int main() {

    //IMPULSE TES
    std::cout << "Impulse test:\n";

    std::vector<cd> impulse = { 1, 0, 0, 0 };
    auto X1 = dft(impulse);

    for (const auto& v : X1)
        std::cout << v << "\n";

    std::cout << "Expected: all (1,0)\n\n";


    //CONSTANT TEST
    std::cout << "Constant signal test:\n";

    std::vector<cd> constant = { 1, 1, 1, 1 };
    auto X2 = dft(constant);

    for (const auto& v : X2)
        std::cout << v << "\n";

    std::cout << "Expected:\n(4,0)\n(0,0)\n(0,0)\n(0,0)\n\n";


    //SINGLE COSINE TEST
    std::cout << "Single-frequency test:\n";

    std::vector<cd> cosine = { 1, 0, -1, 0 };
    auto X3 = dft(cosine);

    for (const auto& v : X3)
        std::cout << v << "\n";

    std::cout << "Expected:\n(0,0)\n(2,0)\n(0,0)\n(2,0)\n\n";


    return 0;
}