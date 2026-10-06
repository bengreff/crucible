// Prototype cost of a 4-D multilinear table lookup (TECHNICAL_PLAN, *Lightweight engine*: the
// budget's "cost per cell update with table lookup"). Not a verification: the table holds smooth
// synthetic values, since only the memory traffic and arithmetic of a lookup are timed.
//
// Layout as planned for the equilibrium table: axes (Z_H, Z_N, e, ln rho), the nv values of one node
// contiguous and ln rho innermost, so a lookup reads 8 runs of 2 nodes. Query streams:
//   random:   uniform over the table (worst case for the caches);
//   coherent: a random walk with small steps, as successive cells of a sweep through a smooth field.
// Usage: crucible_table_lookup_cost [nzh nzn ne nrho nv lookups]
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

namespace {
struct Table {
    std::array<int, 4> n;
    int nv;
    std::vector<double> data;
    double at(int i0, int i1, int i2, int i3, int v) const {
        return data[((((std::size_t)i0 * n[1] + i1) * n[2] + i2) * n[3] + i3) * nv + v];
    }
};
// Multilinear interpolation at unit-cube coordinates x[d] in [0, 1] along each axis.
void lookup(const Table& t, const std::array<double, 4>& x, double* out) {
    std::array<int, 4> i;
    std::array<double, 4> w;
    for (int d = 0; d < 4; ++d) {
        const double s = x[d] * (t.n[d] - 1);
        i[d] = std::min((int)s, t.n[d] - 2);
        w[d] = s - i[d];
    }
    for (int v = 0; v < t.nv; ++v) out[v] = 0;
    for (int c = 0; c < 16; ++c) {
        double weight = 1;
        int j[4];
        for (int d = 0; d < 4; ++d) {
            const int b = (c >> d) & 1;
            j[d] = i[d] + b;
            weight *= b ? w[d] : 1 - w[d];
        }
        const double* p = &t.data[((((std::size_t)j[0] * t.n[1] + j[1]) * t.n[2] + j[2]) * t.n[3] + j[3]) * t.nv];
        for (int v = 0; v < t.nv; ++v) out[v] += weight * p[v];
    }
}
}  // namespace

int main(int argc, char** argv) {
    Table t{{48, 12, 64, 32}, 11, {}};
    long lookups = 2000000;
    if (argc >= 7) {
        for (int d = 0; d < 4; ++d) t.n[d] = std::atoi(argv[1 + d]);
        t.nv = std::atoi(argv[5]);
        lookups = std::atol(argv[6]);
    }
    const std::size_t nodes = (std::size_t)t.n[0] * t.n[1] * t.n[2] * t.n[3];
    t.data.resize(nodes * t.nv);
    for (std::size_t k = 0; k < t.data.size(); ++k) t.data[k] = std::sin(0.001 * k);
    std::printf("table %d x %d x %d x %d nodes, %d values per node, %.1f MB\n", t.n[0], t.n[1], t.n[2], t.n[3], t.nv,
                t.data.size() * 8.0 / 1e6);
    std::mt19937_64 rng(1);
    std::uniform_real_distribution<double> u(0, 1), step(-0.002, 0.002);
    std::vector<double> out(t.nv);
    double sink = 0;
    for (int pattern = 0; pattern < 2; ++pattern) {
        std::vector<std::array<double, 4>> queries(lookups);
        std::array<double, 4> x{0.5, 0.5, 0.5, 0.5};
        for (auto& q : queries) {
            for (int d = 0; d < 4; ++d) {
                x[d] = pattern == 0 ? u(rng) : std::fmin(1, std::fmax(0, x[d] + step(rng)));
            }
            q = x;
        }
        for (int repeat = 0; repeat < 3; ++repeat) {
            const auto start = std::chrono::steady_clock::now();
            for (const auto& q : queries) {
                lookup(t, q, out.data());
                sink += out[0];
            }
            const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            std::printf("%-9s repeat %d: %.1f ns per lookup (%ld lookups)\n", pattern == 0 ? "random" : "coherent",
                        repeat, 1e9 * s / lookups, lookups);
        }
    }
    std::printf("checksum %.6e\n", sink);
}
