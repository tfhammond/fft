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
	std::vector<cd> x = { 1.0, 0.0, 0.0, 0.0 };
	auto X = dft(x);

	for (const auto& v : X) {
		std::cout << v << "\n";
	}
	return 0;
}