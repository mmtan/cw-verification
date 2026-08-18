#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace std;

static void require(bool condition, const string& message) {
    if (!condition) {
        cerr << "VERIFICATION FAILED: " << message << '\n';
        exit(1);
    }
}

struct OrbitData {
    int modulus{};
    vector<vector<int>> orbits;
    vector<int> representatives;
    vector<int> sizes;
    vector<int> orbit_id;
};

static OrbitData multiplier_orbits(int modulus, int multiplier) {
    OrbitData data;
    data.modulus = modulus;
    data.orbit_id.assign(modulus, -1);
    vector<uint8_t> seen(modulus, 0);

    for (int a = 0; a < modulus; ++a) {
        if (seen[a]) continue;
        vector<int> orbit;
        int x = a;
        do {
            orbit.push_back(x);
            seen[x] = 1;
            x = (x * multiplier) % modulus;
        } while (x != a);
        sort(orbit.begin(), orbit.end());
        data.orbits.push_back(orbit);
    }
    sort(data.orbits.begin(), data.orbits.end(),
         [](const vector<int>& left, const vector<int>& right) {
             return left.front() < right.front();
         });

    for (int i = 0; i < static_cast<int>(data.orbits.size()); ++i) {
        data.representatives.push_back(data.orbits[i].front());
        data.sizes.push_back(static_cast<int>(data.orbits[i].size()));
        for (int x : data.orbits[i]) data.orbit_id[x] = i;
    }
    return data;
}

static int correlation(const OrbitData& data,
                       const vector<int>& orbit_values,
                       int shift) {
    int total = 0;
    for (int j = 0; j < data.modulus; ++j) {
        total += orbit_values[data.orbit_id[j]]
               * orbit_values[data.orbit_id[(j + shift) % data.modulus]];
    }
    return total;
}

static vector<int> expand_orbit_vector(const OrbitData& data,
                                       const vector<int>& orbit_values) {
    vector<int> sequence(data.modulus, 0);
    for (int j = 0; j < data.modulus; ++j)
        sequence[j] = orbit_values[data.orbit_id[j]];
    return sequence;
}

static vector<int> decimate(const vector<int>& sequence, int unit) {
    const int modulus = static_cast<int>(sequence.size());
    vector<int> image(modulus, 0);
    for (int j = 0; j < modulus; ++j)
        image[(unit * j) % modulus] = sequence[j];
    return image;
}

static string vector_string(const vector<int>& values) {
    ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << ' ';
        out << values[i];
    }
    return out.str();
}

struct ICWEnumeration {
    uint64_t row_sum_energy_count{};
    vector<vector<int>> solutions;
};

static ICWEnumeration enumerate_multiplier_fixed_icw(const OrbitData& data,
                                                      int coefficient_bound,
                                                      int row_sum,
                                                      int weight) {
    ICWEnumeration answer;
    vector<int> values(data.orbits.size(), 0);
    vector<int> suffix_size(data.orbits.size() + 1, 0);
    for (int i = static_cast<int>(data.orbits.size()) - 1; i >= 0; --i)
        suffix_size[i] = suffix_size[i + 1] + data.sizes[i];

    function<void(int,int,int)> dfs = [&](int index, int sum, int energy) {
        if (energy > weight) return;
        if (index == static_cast<int>(values.size())) {
            if (sum != row_sum || energy != weight) return;
            ++answer.row_sum_energy_count;
            for (int shift = 1; shift < data.modulus; ++shift)
                if (correlation(data, values, shift) != 0) return;
            answer.solutions.push_back(values);
            return;
        }

        const int orbit_size = data.sizes[index];
        const int remaining_size = suffix_size[index + 1];
        for (int value = -coefficient_bound;
             value <= coefficient_bound; ++value) {
            const int next_sum = sum + orbit_size * value;
            if (row_sum < next_sum - coefficient_bound * remaining_size ||
                row_sum > next_sum + coefficient_bound * remaining_size)
                continue;
            const int next_energy = energy + orbit_size * value * value;
            if (next_energy > weight) continue;
            values[index] = value;
            dfs(index + 1, next_sum, next_energy);
        }
    };
    dfs(0, 0, 0);
    return answer;
}

static vector<int> polynomial_multiply_mod3(const vector<int>& left,
                                            const vector<int>& right) {
    vector<int> product(left.size() + right.size() - 1, 0);
    for (size_t i = 0; i < left.size(); ++i)
        for (size_t j = 0; j < right.size(); ++j)
            product[i + j] = (product[i + j] + left[i] * right[j]) % 3;
    for (int& value : product) {
        value %= 3;
        if (value < 0) value += 3;
    }
    return product;
}

static string canonical_shift_and_negation(const array<uint8_t,35>& word) {
    string best;
    bool first = true;
    for (int negate = 0; negate <= 1; ++negate) {
        for (int shift = 0; shift < 35; ++shift) {
            string candidate(35, '0');
            for (int j = 0; j < 35; ++j) {
                int value = word[(j + shift) % 35];
                if (negate && value != 0) value = 3 - value;
                candidate[j] = static_cast<char>('0' + value);
            }
            if (first || candidate < best) {
                best = candidate;
                first = false;
            }
        }
    }
    return best;
}

struct EisensteinInteger {
    int a;
    int b;
}; // a+b*omega, where omega^2+omega+1=0.

static EisensteinInteger sixth_root(int sign, int exponent) {
    static const array<EisensteinInteger,3> roots{{
        {1,0}, {0,1}, {-1,-1}
    }};
    exponent %= 3;
    if (exponent < 0) exponent += 3;
    return {sign * roots[exponent].a, sign * roots[exponent].b};
}

struct SparsePair {
    int shift;
    int left;
    int right;
};

static bool has_eisenstein_unit_lift(const string& ternary_word,
                                     uint64_t& assignments_tested) {
    vector<int> positions;
    vector<int> signs;
    for (int j = 0; j < 35; ++j) {
        const int value = ternary_word[j] - '0';
        if (value != 0) {
            positions.push_back(j);
            signs.push_back(value == 1 ? 1 : -1);
        }
    }
    require(positions.size() == 12,
            "a representative does not have ternary weight 12");

    vector<SparsePair> pairs;
    for (int right = 0; right < 12; ++right) {
        for (int left = 0; left < 12; ++left) {
            if (left == right) continue;
            const int shift = (positions[left] - positions[right] + 35) % 35;
            if (1 <= shift && shift <= 17)
                pairs.push_back({shift, left, right});
        }
    }
    require(pairs.size() == 66,
            "unexpected number of sparse ordered pairs");

    vector<int> phases(12, 0);
    constexpr uint64_t normalized_assignments = 177147; // 3^11.
    for (uint64_t code = 0; code < normalized_assignments; ++code) {
        uint64_t quotient = code;
        phases[0] = 0; // Divide out multiplication by a global power of omega.
        for (int i = 1; i < 12; ++i) {
            phases[i] = static_cast<int>(quotient % 3);
            quotient /= 3;
        }

        array<int,18> coefficient_a{};
        array<int,18> coefficient_b{};
        for (const SparsePair& pair : pairs) {
            EisensteinInteger contribution = sixth_root(
                signs[pair.left] * signs[pair.right],
                phases[pair.left] - phases[pair.right]);
            coefficient_a[pair.shift] += contribution.a;
            coefficient_b[pair.shift] += contribution.b;
        }

        ++assignments_tested;
        bool valid = true;
        for (int shift = 1; shift <= 17; ++shift) {
            if (coefficient_a[shift] != 0 || coefficient_b[shift] != 0) {
                valid = false;
                break;
            }
        }
        if (valid) return true;
    }
    return false;
}

static void verify_common_icw_classification() {
    cout << "=== Multiplier-fixed ICW classification modulo 35 ===\n";
    const OrbitData data = multiplier_orbits(35, 4);
    const vector<int> expected_representatives{0,1,2,3,5,6,7,14,15};
    const vector<int> expected_sizes{1,6,6,6,3,6,2,2,3};
    require(data.representatives == expected_representatives,
            "unexpected multiplier-4 orbit representatives modulo 35");
    require(data.sizes == expected_sizes,
            "unexpected multiplier-4 orbit sizes modulo 35");

    cout << "orbit representatives: " << vector_string(data.representatives) << '\n';
    cout << "orbit sizes: " << vector_string(data.sizes) << '\n';

    const ICWEnumeration bound3 =
        enumerate_multiplier_fixed_icw(data, 3, 6, 36);
    const ICWEnumeration bound4 =
        enumerate_multiplier_fixed_icw(data, 4, 6, 36);

    require(bound3.row_sum_energy_count == 1434,
            "wrong row-sum/energy count for coefficient bound 3");
    require(bound4.row_sum_energy_count == 1600,
            "wrong row-sum/energy count for coefficient bound 4");
    require(bound3.solutions.size() == 2,
            "wrong number of ICW_3(35,36) solutions");
    require(bound4.solutions.size() == 2,
            "wrong number of ICW_4(35,36) solutions");
    require(bound3.solutions == bound4.solutions,
            "the bound-3 and bound-4 solution sets differ");

    const vector<vector<int>> expected_solutions{
        {-3,0,0,0,0,0,0,0,3},
        {-3,0,0,0,3,0,0,0,0}
    };
    require(bound3.solutions == expected_solutions,
            "unexpected invariant ICW solutions");

    const vector<int> first = expand_orbit_vector(data, bound3.solutions[0]);
    const vector<int> second = expand_orbit_vector(data, bound3.solutions[1]);
    require(decimate(first, 3) == second || decimate(second, 3) == first,
            "the two solutions are not equivalent under X -> X^3");

    vector<int> displayed(35, 0);
    displayed[0] = -3;
    displayed[5] = displayed[10] = displayed[20] = 3;
    require(first == displayed || second == displayed,
            "the stated ICW representative was not found");

    cout << "bound 3: row-sum/energy candidates "
         << bound3.row_sum_energy_count << ", full solutions "
         << bound3.solutions.size() << '\n';
    cout << "bound 4: row-sum/energy candidates "
         << bound4.row_sum_energy_count << ", full solutions "
         << bound4.solutions.size() << '\n';
    for (const vector<int>& solution : bound3.solutions)
        cout << "  " << vector_string(solution) << '\n';
    cout << "one equivalence class, represented by\n"
         << "  -3 + 3X^5 + 3X^10 + 3X^20.\n\n";
}

static void verify_cw105_code_obstruction() {
    cout << "=== CW(105,36): ternary-code lifting obstruction ===\n";

    // Low-to-high coefficient lists over F_3.
    const vector<int> x_minus_one{2,1};
    const vector<int> phi5{1,1,1,1,1};
    const vector<int> phi7{1,1,1,1,1,1,1};
    const vector<int> f{
        1,2,2,1,2,1,0,1,2,0,1,0,1
    };
    const vector<int> f_star{
        1,0,1,0,2,1,0,1,2,1,2,2,1
    };
    vector<int> reciprocal(f.rbegin(), f.rend());
    require(reciprocal == f_star,
            "the displayed degree-12 factors are not reciprocal");

    vector<int> factorization{1};
    for (const vector<int>* factor :
         {&x_minus_one, &phi5, &phi7, &f, &f_star})
        factorization = polynomial_multiply_mod3(factorization, *factor);
    vector<int> x35_minus_one(36, 0);
    x35_minus_one[0] = 2;
    x35_minus_one[35] = 1;
    require(factorization == x35_minus_one,
            "the factorization of X^35-1 over F_3 is wrong");

    vector<int> generator{1};
    generator = polynomial_multiply_mod3(generator, x_minus_one);
    generator = polynomial_multiply_mod3(generator, phi5);
    generator = polynomial_multiply_mod3(generator, phi7);
    generator = polynomial_multiply_mod3(generator, f);
    require(generator.size() == 24,
            "the cyclic-code generator should have degree 23");

    map<int,uint64_t> weight_distribution;
    map<string,uint64_t> equivalence_classes;
    constexpr uint64_t messages = 531441; // 3^12.
    for (uint64_t code = 0; code < messages; ++code) {
        array<int,12> message{};
        uint64_t quotient = code;
        for (int i = 0; i < 12; ++i) {
            message[i] = static_cast<int>(quotient % 3);
            quotient /= 3;
        }

        array<uint8_t,35> word{};
        for (int i = 0; i < 12; ++i) {
            if (message[i] == 0) continue;
            for (int j = 0; j < 24; ++j) {
                const int value = word[i + j] + message[i] * generator[j];
                word[i + j] = static_cast<uint8_t>(value % 3);
            }
        }

        int weight = 0;
        for (uint8_t value : word) if (value != 0) ++weight;
        ++weight_distribution[weight];
        if (weight == 12)
            ++equivalence_classes[canonical_shift_and_negation(word)];
    }

    const map<int,uint64_t> expected_distribution{
        {0,1}, {12,420}, {15,2520}, {18,37590},
        {21,158550}, {24,218610}, {27,102620}, {30,11130}
    };
    require(weight_distribution == expected_distribution,
            "unexpected weight distribution of the ternary [35,12] code");
    require(equivalence_classes.size() == 6,
            "the weight-12 words should form six shift/negation classes");
    for (const auto& entry : equivalence_classes)
        require(entry.second == 70,
                "each weight-12 equivalence class should contain 70 words");

    cout << "verified X^35-1 factorization and [35,12] cyclic code\n";
    cout << "weight distribution:\n";
    for (const auto& entry : weight_distribution)
        cout << "  weight " << setw(2) << entry.first << ": "
             << entry.second << '\n';
    cout << "weight-12 words 420; shift/negation classes 6.\n";

    int class_number = 0;
    for (const auto& entry : equivalence_classes) {
        cout << "  class " << ++class_number << ": ";
        for (int j = 0; j < 35; ++j) {
            if (entry.first[j] == '0') continue;
            cout << j << (entry.first[j] == '1' ? '+' : '-') << ' ';
        }
        cout << '\n';
    }

    uint64_t assignments_tested = 0;
    class_number = 0;
    for (const auto& entry : equivalence_classes) {
        const uint64_t before = assignments_tested;
        const bool found =
            has_eisenstein_unit_lift(entry.first, assignments_tested);
        require(!found,
                "an Eisenstein-unit lift satisfying DD^dagger=12 was found");
        cout << "class " << ++class_number << ": tested "
             << assignments_tested - before
             << " normalized phase assignments; solutions 0\n";
    }
    require(assignments_tested == 6ULL * 177147ULL,
            "unexpected total number of phase assignments");
    cout << "total normalized phase assignments: "
         << assignments_tested << '\n';
    cout << "solutions of DD^dagger=12: 0\n\n";
}

static void verify_cw140_local_character_fact() {
    cout << "=== CW(140,36): order-two kernel-character obstruction ===\n";
    int patterns_checked = 0;
    for (int sign : {-1, 1}) {
        const int target_sum = 3 * sign;
        int patterns_for_sign = 0;
        for (int a0 = -1; a0 <= 1; ++a0)
        for (int a1 = -1; a1 <= 1; ++a1)
        for (int a2 = -1; a2 <= 1; ++a2)
        for (int a3 = -1; a3 <= 1; ++a3) {
            if (a0 + a1 + a2 + a3 != target_sum) continue;
            ++patterns_for_sign;
            ++patterns_checked;
            const int alternating_sum = a0 - a1 + a2 - a3;
            require(abs(alternating_sum) == 1,
                    "a four-fibre of sum +/-3 has wrong alternating sum");
            const int energy = a0*a0 + a1*a1 + a2*a2 + a3*a3;
            require(energy == 3,
                    "a four-fibre of sum +/-3 should contain three nonzero entries");
        }
        require(patterns_for_sign == 4,
                "there should be four four-fibre patterns of each sum +/-3");
    }
    require(patterns_checked == 8,
            "unexpected number of local four-fibre patterns");
    cout << "checked all eight four-fibre patterns of sum +/-3;\n"
         << "each has alternating character value +/-1.\n"
         << "The multiplier-fixed ICW_4(35,36) classification above\n"
         << "has only coefficients 0 and +/-3, giving the contradiction.\n\n";
}


struct GaussianInteger {
    int real;
    int imag;
};

static pair<int,int> gaussian_correlation(
        const OrbitData& data,
        const vector<GaussianInteger>& orbit_values,
        int shift) {
    int real = 0;
    int imag = 0;
    for (int j = 0; j < data.modulus; ++j) {
        const GaussianInteger& left =
            orbit_values[data.orbit_id[(j + shift) % data.modulus]];
        const GaussianInteger& right =
            orbit_values[data.orbit_id[j]];
        // (left.real+i left.imag)(right.real-i right.imag).
        real += left.real * right.real + left.imag * right.imag;
        imag += left.imag * right.real - left.real * right.imag;
    }
    return {real, imag};
}

struct GaussianEnumeration {
    uint64_t row_sum_energy_count{};
    vector<uint64_t> filtered_counts;
};

static GaussianEnumeration enumerate_gaussian_orbit_system(
        const OrbitData& data,
        const vector<int>& shifts) {
    GaussianEnumeration answer;
    answer.filtered_counts.assign(shifts.size(), 0);
    vector<GaussianInteger> values(data.orbits.size(), {0,0});
    vector<int> suffix_size(data.orbits.size() + 1, 0);
    for (int i = static_cast<int>(data.orbits.size()) - 1; i >= 0; --i)
        suffix_size[i] = suffix_size[i + 1] + data.sizes[i];

    function<void(int,int,int,int)> dfs =
        [&](int index, int sum_real, int sum_imag, int energy) {
        if (energy > 64) return;
        if (index == static_cast<int>(values.size())) {
            if (sum_real != 8 || sum_imag != 0 || energy != 64) return;
            ++answer.row_sum_energy_count;
            for (size_t k = 0; k < shifts.size(); ++k) {
                if (gaussian_correlation(data, values, shifts[k])
                    != pair<int,int>{0,0})
                    return;
                ++answer.filtered_counts[k];
            }
            return;
        }

        const int orbit_size = data.sizes[index];
        const int remaining_size = suffix_size[index + 1];
        for (int real = -2; real <= 2; ++real) {
            const int next_real = sum_real + orbit_size * real;
            if (8 < next_real - 2 * remaining_size ||
                8 > next_real + 2 * remaining_size)
                continue;
            for (int imag = -2; imag <= 2; ++imag) {
                const int next_imag = sum_imag + orbit_size * imag;
                if (0 < next_imag - 2 * remaining_size ||
                    0 > next_imag + 2 * remaining_size)
                    continue;
                const int next_energy = energy
                    + orbit_size * (real * real + imag * imag);
                if (next_energy > 64) continue;
                values[index] = {real, imag};
                dfs(index + 1, next_real, next_imag, next_energy);
            }
        }
    };
    dfs(0, 0, 0, 0);
    return answer;
}

static void verify_weight64_gaussian_obstructions() {
    cout << "=== Weight 64: Gaussian kernel-character obstructions ===\n";

    struct CaseData {
        int modulus;
        vector<int> representatives;
        vector<int> sizes;
        vector<int> shifts;
        uint64_t initial_count;
        vector<uint64_t> filtered_counts;
    };

    const vector<CaseData> cases{
        {35,
         {0,1,3,5,7,15},
         {1,12,12,3,4,3},
         {1,2,3,4,5},
         1152,
         {4,4,4,4,0}},
        {45,
         {0,1,3,5,7,9,15,21},
         {1,12,4,6,12,4,2,4},
         {1,2,3},
         58188,
         {1242,1242,0}},
        {49,
         {0,1,3,7,21},
         {1,21,21,3,3},
         {1},
         32,
         {0}}
    };

    for (const CaseData& item : cases) {
        const OrbitData data = multiplier_orbits(item.modulus, 2);
        require(data.representatives == item.representatives,
                "unexpected squaring-orbit representatives");
        require(data.sizes == item.sizes,
                "unexpected squaring-orbit sizes");
        const GaussianEnumeration result =
            enumerate_gaussian_orbit_system(data, item.shifts);
        require(result.row_sum_energy_count == item.initial_count,
                "unexpected Gaussian row-sum/energy count");
        require(result.filtered_counts == item.filtered_counts,
                "unexpected Gaussian autocorrelation filtering counts");

        cout << "modulus " << item.modulus << "\n";
        cout << "  orbit representatives: "
             << vector_string(data.representatives) << "\n";
        cout << "  orbit sizes: " << vector_string(data.sizes) << "\n";
        cout << "  row-sum/energy candidates: "
             << result.row_sum_energy_count << "\n";
        for (size_t k = 0; k < item.shifts.size(); ++k)
            cout << "  after Gamma_" << item.shifts[k] << "=0: "
                 << result.filtered_counts[k] << "\n";
        cout << "  full Gaussian solutions: 0\n";
    }
    cout << "The coefficient box [-2,2]+i[-2,2] has no multiplier-fixed\n"
         << "Gaussian sequence of norm 64 for m=35,45,49.\n\n";
}

int main() {
    cout << "Exact verification of the lifting obstructions\n\n";
    require((4 % 35) == ((2 * 2) % 35), "4 is not 2^2 modulo 35");
    int power_of_three = 1;
    for (int i = 0; i < 10; ++i) power_of_three = (3 * power_of_three) % 35;
    require(power_of_three == 4, "4 is not 3^10 modulo 35");

    verify_common_icw_classification();
    verify_cw105_code_obstruction();
    verify_cw140_local_character_fact();
    verify_weight64_gaussian_obstructions();

    cout << "CERTIFICATE VERIFIED:\n"
         << "  the exact finite obstructions used for CW(105,36), CW(140,36),\n"
         << "  CW(140,64), CW(180,64), and CW(196,64) have been reproduced\n"
         << "  successfully.\n";
    return 0;
}
