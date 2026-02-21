#include <vector>
#include <complex>
#include <cmath>
#include <numbers>
#include <iostream>

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