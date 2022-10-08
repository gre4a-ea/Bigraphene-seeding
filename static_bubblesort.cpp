#include <iostream>

int main() {
    const auto size = 1000U; 
    
    int a[size];
    
    std::size_t n; //unsigned int
    std::cin >> n;
    
    if (n > size) {
		std::cerr << "Too many elements.\n";
		return EXIT_FAILURE;
	}
	
	for (std::size_t i = 0; i < n; ++i) {
		std::cin >> a[i];
	}
	
    for (std::size_t i = 0; i < n - 1; ++i) {
		for (std::size_t j = i + 1; j < n; ++j) {
			if (a[i] > a[j]) {
				std::swap(a[i], a[j]);
			}
		}
	}
    
    for (std::size_t i = 0; i < n; ++i) {
		std::cout << a[i] << ' ';
	}
    
    return EXIT_SUCCESS;
    
}
