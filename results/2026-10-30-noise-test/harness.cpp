// Three-arm integer harness around the repo's simulator (included unchanged). See PREREG.md.
#define main repo_main
#include "asolaria_taylor_green_regularity.cpp"
#undef main
#include <cstdint>
using i64 = long long;

static int sgn_s(int i) { return i < 12 ? 1 : -1; }
static int sgn_c(int i) { return (i < 6 || i >= 18) ? 1 : -1; }

struct Row { int t; i64 e2, kj, wj, a; int occ_bad; };

static void step(vector<u8>& state, int t, int L, const vector<Collision>& c, const array<Velocity, K>& v) {
    int V = L * L * L;
    const int count = 32;
    const size_t total_c = c.size();
    for (int cell = 0; cell < V; ++cell) {
        int base = cell * K;
        for (int j = 0; j < count; ++j) {
            size_t which = (static_cast<size_t>(t) * count + static_cast<size_t>(cell) * 17 + static_cast<size_t>(j) * 37) % total_c;
            const auto& col = c[which];
            int code = encode(state[base + col.a], state[base + col.b], state[base + col.c], state[base + col.d]);
            auto dest = decode(col.lut[code]);
            state[base + col.a] = static_cast<u8>(dest[0]);
            state[base + col.b] = static_cast<u8>(dest[1]);
            state[base + col.c] = static_cast<u8>(dest[2]);
            state[base + col.d] = static_cast<u8>(dest[3]);
        }
    }
    vector<u8> next(V * K);
    for (int x = 0; x < L; ++x) for (int y = 0; y < L; ++y) for (int z = 0; z < L; ++z)
        for (int ki = 0; ki < K; ++ki) {
            int nx = (x + v[ki].x + L) % L, ny = (y + v[ki].y + L) % L, nz = (z + v[ki].z + L) % L;
            next[((nx * L + ny) * L + nz) * K + ki] = state[((x * L + y) * L + z) * K + ki];
        }
    state.swap(next);
}

static Row measure(const vector<u8>& s, int t, int L, const array<Velocity, K>& v) {
    int V = L * L * L;
    vector<int> jx(V), jy(V), jz(V);
    Row r{t, 0, 0, 0, 0, 0};
    for (int cell = 0; cell < V; ++cell) {
        int a = 0, b = 0, d = 0;
        for (int ki = 0; ki < K; ++ki) {
            int n = s[cell * K + ki];
            if (n > 2) r.occ_bad++;
            a += n * v[ki].x; b += n * v[ki].y; d += n * v[ki].z; r.e2 += n * v[ki].e;
        }
        jx[cell] = a; jy[cell] = b; jz[cell] = d;
        r.kj += (i64)a * a + (i64)b * b + (i64)d * d;
    }
    auto id = [&](int x, int y, int z) { x = (x % L + L) % L; y = (y % L + L) % L; z = (z % L + L) % L; return (x * L + y) * L + z; };
    for (int x = 0; x < L; ++x) for (int y = 0; y < L; ++y) for (int z = 0; z < L; ++z) {
        i64 wx = (jz[id(x, y + 1, z)] - jz[id(x, y - 1, z)]) - (jy[id(x, y, z + 1)] - jy[id(x, y, z - 1)]);
        i64 wy = (jx[id(x, y, z + 1)] - jx[id(x, y, z - 1)]) - (jz[id(x + 1, y, z)] - jz[id(x - 1, y, z)]);
        i64 wz = (jy[id(x + 1, y, z)] - jy[id(x - 1, y, z)]) - (jx[id(x, y + 1, z)] - jx[id(x, y - 1, z)]);
        r.wj += wx * wx + wy * wy + wz * wz;
        int cell = id(x, y, z);
        r.a += (i64)jx[cell] * sgn_s(x) * sgn_c(y) * sgn_c(z) - (i64)jy[cell] * sgn_c(x) * sgn_s(y) * sgn_c(z);
    }
    return r;
}

int main() {
    const int L = 24, V = L * L * L, T = 480;
    auto v = make_velocities();
    auto c = make_collisions(v);
    vector<u8> tgv(V * K), nul(V * K);
    std::streambuf* old = std::cout.rdbuf(nullptr);  // silence the repo initialiser's prints
    init_taylor_green_exact(L, 0.20, tgv, v);
    init_taylor_green_exact(L, 0.0, nul, v);
    std::cout.rdbuf(old);
    // SCR: Fisher-Yates over whole cell states, integer LCG
    vector<int> perm(V);
    for (int i = 0; i < V; ++i) perm[i] = i;
    uint64_t x = 20261030;
    for (int i = V - 1; i > 0; --i) {
        x = (1103515245ULL * x + 12345ULL) % 2147483648ULL;
        int j = (int)(x % (uint64_t)(i + 1));
        std::swap(perm[i], perm[j]);
    }
    vector<u8> scr(V * K);
    for (int i = 0; i < V; ++i) for (int ki = 0; ki < K; ++ki) scr[i * K + ki] = tgv[perm[i] * K + ki];
    const char* names[3] = {"TGV", "SCR", "NULL"};
    vector<u8>* st[3] = {&tgv, &scr, &nul};
    vector<vector<Row>> rows(3);
    for (int a = 0; a < 3; ++a)
        for (int t = 0; t <= T; ++t) {
            rows[a].push_back(measure(*st[a], t, L, v));
            if (t < T) step(*st[a], t, L, c, v);
        }
    for (int a = 0; a < 3; ++a)
        for (auto& r : rows[a])
            printf("TGROW|arm=%s|t=%d|E2=%lld|KJ=%lld|WJ=%lld|A=%lld|occ_bad=%d|json=0\n", names[a], r.t, r.e2, r.kj, r.wj, r.a, r.occ_bad);
    // predictions
    bool p1 = true, p2 = true;
    for (int a = 0; a < 3; ++a) for (auto& r : rows[a]) p1 &= r.occ_bad == 0 && r.e2 == rows[a][0].e2;
    for (auto& r : rows[2]) p2 &= r.kj == 0 && r.wj == 0 && r.a == 0;
    i64 a0 = rows[0][0].a, as0 = rows[1][0].a;
    bool p3 = a0 > 0 && (as0 < 0 ? -as0 : as0) * 10 <= a0;
    int within = 0, above = 0, n = 0;
    for (int t = 5; t <= T; ++t) {
        i64 wt = rows[0][t].wj, ws = rows[1][t].wj, d = wt > ws ? wt - ws : ws - wt;
        n++; if (d * 10 <= ws) within++; if (wt * 10 > ws * 11) above++;
    }
    const char* p4 = (within * 10 >= n * 9) ? "NOISE" : (above * 10 >= n * 9) ? "VORTEX" : "MIXED";
    int half = -1;
    for (int t = 0; t <= T; ++t) if (rows[0][t].a * 2 <= a0) { half = t; break; }
    printf("TGSUM|ticks_5_480=%d|within_10pct=%d|tgv_above_110pct=%d|P4=%s|json=0\n", n, within, above, p4);
    printf("TGSUM|A_TGV0=%lld|A_SCR0=%lld|A_TGV_half_tick=%d|WJ_TGV0=%lld|WJ_TGV1=%lld|WJ_SCR0=%lld|WJ_SCR1=%lld|json=0\n",
           a0, as0, half, rows[0][0].wj, rows[0][1].wj, rows[1][0].wj, rows[1][1].wj);
    printf("TGRESULT|P1=%s|P2=%s|P3=%s|P4=%s|json=0\n", p1 ? "PASS" : "FAIL", p2 ? "PASS" : "FAIL", p3 ? "PASS" : "FAIL", p4);
}
