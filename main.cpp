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

// c) Тест на самую длинную последовательность единиц в блоке
double longest_run_of_ones_test(const std::bitset<128>& sequence) {
    int k = 5;
    int n = sequence.size();
    int m = 8;

    std::vector<int> counts(6, 0);
    int consecutive_ones = 0;

    for (size_t i = 0; i < sequence.size(); ++i) {
        if (sequence[i]) {
            consecutive_ones++;
            if (i == sequence.size() - 1 && consecutive_ones >= 1 && consecutive_ones <= 5) {
                counts[consecutive_ones - 1]++;
            }
        } else {
            if (consecutive_ones >= 1 && consecutive_ones <= 5) {
                counts[consecutive_ones - 1]++;
            }
            consecutive_ones = 0;
        }
    }

    std::vector<double> pi = {0.2148, 0.3672, 0.2305, 0.1250, 0.0625};
    double chi_square = 0.0;

    for (int i = 0; i < k; ++i) {
        chi_square += std::pow((counts[i] - n * pi[i] / m), 2) / (n * pi[i] / m);
    }

    double p_value = std::exp(-chi_square / 2);

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
