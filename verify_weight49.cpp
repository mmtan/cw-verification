#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using std::array;
using std::cerr;
using std::cout;
using std::endl;
using std::int64_t;
using std::string;
using std::vector;

struct OrbitData {
    int modulus;
    int multiplier;
    vector<vector<int>> orbits;
    vector<int> orbit_id;
    vector<int> reps;
    vector<int> sizes;
};

static OrbitData make_orbits(int modulus, int multiplier) {
    OrbitData data;
    data.modulus = modulus;
    data.multiplier = multiplier;
    data.orbit_id.assign(modulus, -1);
    vector<char> seen(modulus, 0);
    for (int start = 0; start < modulus; ++start) {
        if (seen[start]) continue;
        vector<int> orbit;
        int x = start;
        while (!seen[x]) {
            seen[x] = 1;
            orbit.push_back(x);
            x = (x * multiplier) % modulus;
        }
        int rep = *std::min_element(orbit.begin(), orbit.end());
        data.orbits.push_back(std::move(orbit));
        data.reps.push_back(rep);
    }
    // The scan by start already orders orbits by their minimum element, but sort
    // explicitly so that the output and certificate are deterministic.
    vector<int> order(data.orbits.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return data.reps[a] < data.reps[b];
    });
    vector<vector<int>> sorted_orbits;
    vector<int> sorted_reps;
    for (int idx : order) {
        sorted_orbits.push_back(data.orbits[idx]);
        sorted_reps.push_back(data.reps[idx]);
    }
    data.orbits = std::move(sorted_orbits);
    data.reps = std::move(sorted_reps);
    data.sizes.clear();
    std::fill(data.orbit_id.begin(), data.orbit_id.end(), -1);
    for (int ell = 0; ell < static_cast<int>(data.orbits.size()); ++ell) {
        data.sizes.push_back(static_cast<int>(data.orbits[ell].size()));
        for (int j : data.orbits[ell]) data.orbit_id[j] = ell;
    }
    return data;
}

static int correlation(const OrbitData& data, const vector<int8_t>& x, int shift) {
    int total = 0;
    for (int j = 0; j < data.modulus; ++j) {
        total += static_cast<int>(x[data.orbit_id[j]]) *
                 static_cast<int>(x[data.orbit_id[(j + shift) % data.modulus]]);
    }
    return total;
}

static void print_vector(const vector<int8_t>& x) {
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (i) cout << ' ';
        cout << static_cast<int>(x[i]);
    }
}

static bool equal_vec(const vector<int>& a, const vector<int>& b) {
    return a == b;
}

static void require(bool cond, const string& message) {
    if (!cond) {
        cerr << "CERTIFICATE CHECK FAILED: " << message << '\n';
        std::exit(1);
    }
}

static void verify_116() {
    cout << "=== CW(116,49) orbit certificate ===\n";
    const OrbitData data = make_orbits(116, 7);
    const vector<int> expected_reps = {0,1,2,3,4,5,6,8,10,15,16,29,30,32,58};
    const vector<int> expected_sizes = {1,14,7,14,7,14,7,7,7,14,7,2,7,7,1};
    require(equal_vec(data.reps, expected_reps), "unexpected orbit representatives modulo 116");
    require(equal_vec(data.sizes, expected_sizes), "unexpected orbit sizes modulo 116");

    cout << "orbit representatives:";
    for (int r : data.reps) cout << ' ' << r;
    cout << "\norbit sizes:";
    for (int w : data.sizes) cout << ' ' << w;
    cout << '\n';

    std::array<int64_t, 5> count{}; // sum/energy, then C1,...,C4
    vector<vector<int8_t>> before_c4;
    vector<int> c4_values;
    vector<int8_t> x(data.orbits.size(), 0);

    std::function<void(int,int,int)> dfs = [&](int idx, int sum, int energy) {
        if (energy > 49) return;
        if (idx == static_cast<int>(x.size())) {
            if (sum != 7 || energy != 49) return;
            ++count[0];
            if (correlation(data, x, 1) != 0) return;
            ++count[1];
            if (correlation(data, x, 2) != 0) return;
            ++count[2];
            if (correlation(data, x, 3) != 0) return;
            ++count[3];
            before_c4.push_back(x);
            int c4 = correlation(data, x, 4);
            c4_values.push_back(c4);
            if (c4 == 0) ++count[4];
            return;
        }
        int remaining_weight = 0;
        for (int j = idx + 1; j < static_cast<int>(x.size()); ++j) remaining_weight += data.sizes[j];
        const int w = data.sizes[idx];
        for (int value = -1; value <= 1; ++value) {
            int sum2 = sum + w * value;
            // Remaining entries lie in {-1,0,1}.
            if (7 < sum2 - remaining_weight || 7 > sum2 + remaining_weight) continue;
            int energy2 = energy + w * value * value;
            if (energy2 > 49) continue;
            x[idx] = static_cast<int8_t>(value);
            dfs(idx + 1, sum2, energy2);
        }
    };
    dfs(0,0,0);

    const std::array<int64_t,5> expected = {6088,624,40,16,0};
    require(count == expected, "unexpected filtering counts for modulus 116");
    require(before_c4.size() == 16, "unexpected survivor count before C4");
    for (int value : c4_values) require(value == -1 || value == -3, "unexpected C4 value");

    cout << "row sum and energy: " << count[0] << '\n';
    cout << "after C1=0:       " << count[1] << '\n';
    cout << "after C2=0:       " << count[2] << '\n';
    cout << "after C3=0:       " << count[3] << '\n';
    cout << "after C4=0:       " << count[4] << '\n';
    cout << "survivors before C4, followed by C4:\n";
    for (std::size_t i = 0; i < before_c4.size(); ++i) {
        print_vector(before_c4[i]);
        cout << " | " << c4_values[i] << '\n';
    }
    cout << "certificate: every C1=C2=C3=0 survivor has C4 in {-1,-3}.\n\n";
}

static void verify_120_folded() {
    cout << "=== Folded obstruction for CW(20q,49), q<=6 ===\n";
    const OrbitData data = make_orbits(20, 7);
    const vector<int> expected_reps = {0,1,2,4,5,10,11};
    const vector<int> expected_sizes = {1,4,4,4,2,1,4};
    require(equal_vec(data.reps, expected_reps), "unexpected orbit representatives modulo 20");
    require(equal_vec(data.sizes, expected_sizes), "unexpected orbit sizes modulo 20");
    // Exhaustively check the strongest coefficient box [-6,6].  A solution for
    // any q<=6 would occur in this box.
    vector<int8_t> x(data.orbits.size(), 0);
    int64_t sum_energy = 0;
    int64_t full = 0;
    std::function<void(int,int,int)> dfs = [&](int idx, int sum, int energy) {
        if (energy > 49) return;
        if (idx == static_cast<int>(x.size())) {
            if (sum != 7 || energy != 49) return;
            ++sum_energy;
            bool ok = true;
            for (int t = 1; t < 20; ++t) {
                if (correlation(data, x, t) != 0) { ok = false; break; }
            }
            if (ok) ++full;
            return;
        }
        int remaining_weight = 0;
        for (int j = idx + 1; j < static_cast<int>(x.size()); ++j) remaining_weight += data.sizes[j];
        const int w = data.sizes[idx];
        for (int value = -6; value <= 6; ++value) {
            int sum2 = sum + w * value;
            if (7 < sum2 - 6 * remaining_weight || 7 > sum2 + 6 * remaining_weight) continue;
            int energy2 = energy + w * value * value;
            if (energy2 > 49) continue;
            x[idx] = static_cast<int8_t>(value);
            dfs(idx + 1, sum2, energy2);
        }
    };
    dfs(0,0,0);
    require(full == 0, "unexpected folded solution modulo 20");
    cout << "orbit representatives:";
    for (int r : data.reps) cout << ' ' << r;
    cout << "\norbit sizes:";
    for (int w : data.sizes) cout << ' ' << w;
    cout << "\nrow-sum/energy candidates in [-6,6]^7: " << sum_energy;
    cout << "\nfull folded solutions: " << full << "\n\n";
}

struct FeasibilityDP {
    // feasible[idx][energy][sum + offset] for suffix idx..end
    int offset;
    int max_sum;
    vector<vector<vector<uint8_t>>> feasible;
};

static FeasibilityDP build_feasibility(const vector<int>& weights, int bound, int target_energy) {
    const int total_weight = std::accumulate(weights.begin(), weights.end(), 0);
    const int max_sum = bound * total_weight;
    const int offset = max_sum;
    const int t = static_cast<int>(weights.size());
    FeasibilityDP dp;
    dp.offset = offset;
    dp.max_sum = max_sum;
    dp.feasible.assign(t + 1,
        vector<vector<uint8_t>>(target_energy + 1,
            vector<uint8_t>(2 * max_sum + 1, 0)));
    dp.feasible[t][0][offset] = 1;
    for (int idx = t - 1; idx >= 0; --idx) {
        int w = weights[idx];
        for (int e = 0; e <= target_energy; ++e) {
            for (int s = -max_sum; s <= max_sum; ++s) {
                bool ok = false;
                for (int value = -bound; value <= bound && !ok; ++value) {
                    int e_rem = e - w * value * value;
                    int s_rem = s - w * value;
                    if (e_rem < 0 || std::abs(s_rem) > max_sum) continue;
                    if (dp.feasible[idx + 1][e_rem][s_rem + offset]) ok = true;
                }
                dp.feasible[idx][e][s + offset] = static_cast<uint8_t>(ok);
            }
        }
    }
    return dp;
}

static void verify_192() {
    cout << "=== CW(192,49) folded orbit certificate ===\n";
    const OrbitData data = make_orbits(64, 7);
    const vector<int> expected_reps = {0,1,2,3,4,6,8,9,11,12,16,18,22,24,32,36,44};
    const vector<int> expected_sizes = {1,8,4,8,2,4,2,8,8,2,2,4,4,2,1,2,2};
    require(equal_vec(data.reps, expected_reps), "unexpected orbit representatives modulo 64");
    require(equal_vec(data.sizes, expected_sizes), "unexpected orbit sizes modulo 64");

    cout << "orbit representatives:";
    for (int r : data.reps) cout << ' ' << r;
    cout << "\norbit sizes:";
    for (int w : data.sizes) cout << ' ' << w;
    cout << '\n';

    const vector<int> shifts = {1,2,3,4,5,6,7,8,11,12,16};
    vector<int64_t> counts(shifts.size() + 1, 0); // initial + after each shift
    vector<vector<int8_t>> survivors;
    vector<int> c24_values;
    vector<int8_t> x(data.orbits.size(), 0);
    const FeasibilityDP dp = build_feasibility(data.sizes, 3, 49);

    std::function<void(int,int,int)> dfs = [&](int idx, int sum, int energy) {
        if (idx == static_cast<int>(x.size())) {
            if (sum != 7 || energy != 49) return;
            ++counts[0];
            for (std::size_t k = 0; k < shifts.size(); ++k) {
                if (correlation(data, x, shifts[k]) != 0) return;
                ++counts[k + 1];
            }
            survivors.push_back(x);
            c24_values.push_back(correlation(data, x, 24));
            return;
        }
        const int w = data.sizes[idx];
        for (int value = -3; value <= 3; ++value) {
            int sum2 = sum + w * value;
            int energy2 = energy + w * value * value;
            if (energy2 > 49) continue;
            int need_sum = 7 - sum2;
            int need_energy = 49 - energy2;
            if (std::abs(need_sum) > dp.max_sum) continue;
            if (!dp.feasible[idx + 1][need_energy][need_sum + dp.offset]) continue;
            x[idx] = static_cast<int8_t>(value);
            dfs(idx + 1, sum2, energy2);
        }
    };
    dfs(0,0,0);

    const vector<int64_t> expected = {
        22880810,
        5952866,
        455924,
        413640,
        14670,
        14670,
        7836,
        7836,
        338,
        338,
        22,
        4
    };
    require(counts == expected, "unexpected filtering counts for folded modulus 64");
    require(survivors.size() == 4, "unexpected survivor count before C24");
    for (int value : c24_values) require(value == 12, "unexpected C24 value");

    cout << "row sum and energy: " << counts[0] << '\n';
    for (std::size_t k = 0; k < shifts.size(); ++k) {
        cout << "after C" << shifts[k] << "=0:";
        int width = shifts[k] < 10 ? 12 : 11;
        cout << std::setw(width) << counts[k + 1] << '\n';
    }
    cout << "survivors before C24, followed by C24:\n";
    for (std::size_t i = 0; i < survivors.size(); ++i) {
        print_vector(survivors[i]);
        cout << " | " << c24_values[i] << '\n';
    }
    cout << "certificate: every listed survivor has C24=12.\n\n";
}

int main() {
    auto start = std::chrono::steady_clock::now();
    verify_120_folded();
    verify_116();
    verify_192();
    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    cout << "All exact certificate checks passed.\n";
    cout << std::fixed << std::setprecision(3)
         << "elapsed seconds: " << elapsed.count() << '\n';
    return 0;
}
