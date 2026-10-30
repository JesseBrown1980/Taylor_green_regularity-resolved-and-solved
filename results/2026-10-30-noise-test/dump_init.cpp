// Exports the subject's start states and gate tables (it uses the repo's own float initialiser).
#define main repo_main
#include "asolaria_taylor_green_regularity.cpp"
#undef main
int main() {
    const int L = 24, V = L * L * L;
    auto v = make_velocities();
    auto c = make_collisions(v);
    vector<u8> tgv(V * K), nul(V * K);
    std::streambuf* old = std::cout.rdbuf(nullptr);
    init_taylor_green_exact(L, 0.20, tgv, v);
    init_taylor_green_exact(L, 0.0, nul, v);
    std::cout.rdbuf(old);
    std::ofstream("INIT_TGV.trits", std::ios::binary).write((const char*)tgv.data(), tgv.size());
    std::ofstream("INIT_NULL.trits", std::ios::binary).write((const char*)nul.data(), nul.size());
    std::ofstream g("GATES.lut", std::ios::binary);
    for (auto& col : c) { u8 h[4] = {(u8)col.a, (u8)col.b, (u8)col.c, (u8)col.d}; g.write((const char*)h, 4); g.write((const char*)col.lut.data(), 81); }
    printf("gates=%zu\n", c.size());
}
