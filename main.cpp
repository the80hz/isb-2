#include <iostream>
#include <bitset>
#include <random>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

std::vector<std::bitset<128>> generate_random_sequences(size_t count) {
    std::random_device rd;
    std::mt19937_64 generator(rd());
    std::uniform_int_distribution<uint64_t> distribution(0, UINT64_MAX);

    std::vector<std::bitset<128>> sequences;
    sequences.reserve(count);

    for (size_t i = 0; i < count; ++i) {
        uint64_t first_half = distribution(generator);
        uint64_t second_half = distribution(generator);
        std::bitset<64> first_half_bitset(first_half);
        std::bitset<64> second_half_bitset(second_half);
        std::string full_bit_string = first_half_bitset.to_string() + second_half_bitset.to_string();
        std::bitset<128> sequence(full_bit_string);
        sequences.push_back(sequence);
    }

    return sequences;
}


// a) Частотный побитовый тест
double monobit_frequency_test(const std::bitset<128>& sequence) {
    int sum = 0;

    for (size_t i = 0; i < sequence.size(); ++i) {
        sum += sequence[i] ? 1 : -1;
    }

    double s_obs = std::abs(sum) / std::sqrt(sequence.size());
    double p_value = std::erfc(s_obs / std::sqrt(2));

    return p_value;
}

// b) Тест на одинаковые подряд идущие биты
double runs_test(const std::bitset<128>& sequence) {
    size_t runs_count = 0;
    size_t ones_count = 0;

    for (size_t i = 0; i < sequence.size(); ++i) {
        if (sequence[i]) {
            ++ones_count;
            if (i == 0 || !sequence[i - 1]) {
                ++runs_count;
            }
        } else {
            if (i == 0 || sequence[i - 1]) {
                ++runs_count;
            }
        }
    }

    double pi = static_cast<double>(ones_count) / sequence.size();
    double tau = 2.0 / std::sqrt(sequence.size());
    if (std::abs(pi - 0.5) >= tau) {
        return 0.0;
    }

    double test_statistic = (runs_count - (2.0 * sequence.size() * pi * (1.0 - pi))) / std::sqrt(2.0 * sequence.size() * pi * (1.0 - pi));
    double p_value = std::erfc(std::abs(test_statistic) / std::sqrt(2.0));
    return p_value;
}

double incomplete_gamma_integral(double s, double x) {
    auto integrand = [s](double t) {
        return std::pow(t, s - 1) * std::exp(-t);
    };

    const int n = 50;
    const double a = 0;
    const double b = x;
    const double h = (b - a) / (2 * n);

    double result = 0.0;
    double x_i = a + h;

    for (int i = 1; i <= n; ++i) {
        result += integrand(x_i - h / 2);
        result += 2 * integrand(x_i);
        x_i += h;
    }

    result += integrand(b - h / 2);
    result *= h / 3;

    return result;
}

double chi_squared_cdf(double x, double k) {
    return incomplete_gamma_integral(k / 2.0, x / 2.0) / std::tgamma(k / 2.0);
}

// c) Тест на самую длинную последовательность единиц в блоке
double longest_run_of_ones_test(const std::bitset<128>& sequence) {
    std::vector<size_t> frequencies(6, 0);
    size_t current_run = 0;

    for (size_t i = 0; i < sequence.size(); ++i) {
        if (sequence[i]) {
            ++current_run;
        } else {
            if (current_run > 0) {
                if (current_run <= 6) {
                    ++frequencies[current_run - 1];
                } else {
                    ++frequencies[5];
                }
                current_run = 0;
            }
        }
    }

    if (current_run > 0) {
        if (current_run <= 6) {
            ++frequencies[current_run - 1];
        } else {
            ++frequencies[5];
        }
    }

    double chi_squared = 0.0;
    std::vector<double> probabilities = {0.2148, 0.3672, 0.2305, 0.1250, 0.0463, 0.0162};
    for (size_t i = 0; i < 6; ++i) {
        double diff = static_cast<double>(frequencies[i]) - probabilities[i] * sequence.size();
        chi_squared += (diff * diff) / (probabilities[i] * sequence.size());
    }

    double p_value = 1.0 - chi_squared_cdf(chi_squared, 5);
    return p_value;
}

int main() {
    size_t sequence_count = 4;
    std::vector<std::bitset<128>> sequences = generate_random_sequences(sequence_count);

    for (size_t i = 0; i < sequences.size(); ++i) {
        const auto& sequence = sequences[i];
        std::cout << "Sequence " << (i + 1) << ":\n";
        std::cout << sequence << std::endl;

        double p_value_a = monobit_frequency_test(sequence);
        double p_value_b = runs_test(sequence);
        double p_value_c = longest_run_of_ones_test(sequence);

        std::cout << "a) Test for frequency bitwise: p_value = " << p_value_a << std::endl;
        std::cout << "b) Test for identical consecutive bits: p_value = " << p_value_b << std::endl;
        std::cout << "c) Test for the longest sequence of units in a block: p_value = " << p_value_c << std::endl;

        std::cout << std::endl;
    }

    return 0;
}
