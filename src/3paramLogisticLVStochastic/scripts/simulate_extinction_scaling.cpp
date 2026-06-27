#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);
    return 0;
}();

struct Env
{
    double r;
    double K;
    double gamma;
};

struct State
{
    double x;
    double logy;
};

constexpr Env ENV_A{1.0, 1.0, 1.0};
constexpr Env ENV_B{4.0, 2.0, 2.2};

constexpr double DEFAULT_T_MAX = 12000.0;
constexpr long long DEFAULT_N_TRAJ = 20000;
constexpr double DEFAULT_DT = 0.01;
constexpr double LOG_Y_MIN = -700.0;
constexpr double LOG_Y_MAX = 700.0;
constexpr double X_MIN = 1e-14;

// Barrier sweep: extinction thresholds y_ext, geometric from Y_HI down to Y_LO.
// One trajectory is scored against ALL of them at once via its running minimum,
// so the whole tau-vs-barrier curve costs a single ensemble pass.
constexpr double Y_HI = 1e-1;
constexpr double Y_LO = 1e-5;
constexpr int N_LEVELS = 17;

inline double safe_exp(double value)
{
    if (value <= LOG_Y_MIN)
        return 0.0;
    if (value >= LOG_Y_MAX)
        return std::exp(LOG_Y_MAX);
    return std::exp(value);
}

inline double clamp_positive(double value, double floor_value)
{
    return (value > floor_value) ? value : floor_value;
}

inline void rhs(const Env &env, const State &state, double &dx, double &dlogy)
{
    double y = safe_exp(state.logy);
    dx = env.r * state.x * (1.0 - state.x / env.K) - state.x * y;
    dlogy = state.x - env.gamma;
}

inline void rk4_step(const Env &env, State &state, double dt)
{
    double k1x = 0.0, k1logy = 0.0;
    rhs(env, state, k1x, k1logy);
    State s2{state.x + 0.5 * dt * k1x, state.logy + 0.5 * dt * k1logy};
    double k2x = 0.0, k2logy = 0.0;
    rhs(env, s2, k2x, k2logy);
    State s3{state.x + 0.5 * dt * k2x, state.logy + 0.5 * dt * k2logy};
    double k3x = 0.0, k3logy = 0.0;
    rhs(env, s3, k3x, k3logy);
    State s4{state.x + dt * k3x, state.logy + dt * k3logy};
    double k4x = 0.0, k4logy = 0.0;
    rhs(env, s4, k4x, k4logy);
    state.x += (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    state.logy += (dt / 6.0) * (k1logy + 2.0 * k2logy + 2.0 * k3logy + k4logy);
    state.x = clamp_positive(state.x, X_MIN);
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " alpha [t_max] [n_traj] [dt] [seed]\n";
        return 1;
    }

    const double alpha = std::stod(argv[1]);
    const double t_max = (argc > 2) ? std::stod(argv[2]) : DEFAULT_T_MAX;
    const long long n_traj = (argc > 3) ? std::stoll(argv[3]) : DEFAULT_N_TRAJ;
    const double dt = (argc > 4) ? std::stod(argv[4]) : DEFAULT_DT;
    const std::uint64_t base_seed =
        (argc > 5) ? static_cast<std::uint64_t>(std::stoull(argv[5])) : 0x9E3779B97F4A7C15ULL;

    if (alpha <= 0.0 || t_max <= 0.0 || n_traj <= 0 || dt <= 0.0)
    {
        std::cerr << "Invalid parameters.\n";
        return 1;
    }

    const double x0 = 0.5 * (ENV_A.K + ENV_B.K);
    const double y0 = 0.5 * (ENV_A.K + ENV_B.K);

    // Barrier levels in log-space, sorted from shallowest (highest y_ext) to
    // deepest (lowest y_ext), i.e. decreasing log-threshold.
    std::vector<double> y_ext(N_LEVELS);
    std::vector<double> log_level(N_LEVELS);
    for (int k = 0; k < N_LEVELS; ++k)
    {
        double frac = static_cast<double>(k) / (N_LEVELS - 1);
        y_ext[k] = Y_HI * std::pow(Y_LO / Y_HI, frac); // geometric, decreasing
        log_level[k] = std::log(y_ext[k]);
    }

    unsigned int hw = std::thread::hardware_concurrency();
    unsigned int n_threads = (hw > 3) ? hw - 2 : 1;
    if (static_cast<long long>(n_threads) > n_traj)
        n_threads = static_cast<unsigned int>(n_traj);

    // Per-thread accumulators, combined after the join (no shared writes).
    std::vector<std::vector<long long>> n_died(n_threads, std::vector<long long>(N_LEVELS, 0));
    std::vector<std::vector<double>> exposure(n_threads, std::vector<double>(N_LEVELS, 0.0));

    auto worker = [&](unsigned int tid, long long lo, long long hi)
    {
        std::vector<long long> &nd = n_died[tid];
        std::vector<double> &ex = exposure[tid];
        for (long long i = lo; i < hi; ++i)
        {
            std::uint64_t seed = base_seed + 0x100000001B3ULL * static_cast<std::uint64_t>(i);
            std::mt19937_64 rng(seed);
            std::exponential_distribution<double> dwell(alpha);
            std::bernoulli_distribution coin(0.5);

            State state{x0, std::log(y0)};
            const Env *env = coin(rng) ? &ENV_A : &ENV_B;

            // First-passage time to each level for this trajectory (-1 = never).
            // Levels are crossed in order of increasing depth as logy descends.
            int next_k = 0;
            double cross_time[N_LEVELS];
            for (int k = 0; k < N_LEVELS; ++k)
                cross_time[k] = -1.0;

            double t = 0.0;
            while (t < t_max && next_k < N_LEVELS)
            {
                double dwell_end = std::min(t + dwell(rng), t_max);
                while (t < dwell_end)
                {
                    double step = std::min(dt, dwell_end - t);
                    rk4_step(*env, state, step);
                    t += step;
                    while (next_k < N_LEVELS && state.logy <= log_level[next_k])
                    {
                        cross_time[next_k] = t;
                        ++next_k;
                    }
                    if (next_k >= N_LEVELS)
                        break;
                }
                env = (env == &ENV_A) ? &ENV_B : &ENV_A;
            }

            // Tally: a level that was crossed contributes its crossing time and
            // one death; an uncrossed level is censored at t_max.
            for (int k = 0; k < N_LEVELS; ++k)
            {
                if (cross_time[k] >= 0.0)
                {
                    ++nd[k];
                    ex[k] += cross_time[k];
                }
                else
                {
                    ex[k] += t_max;
                }
            }
        }
    };

    std::vector<std::thread> pool;
    long long chunk = (n_traj + n_threads - 1) / static_cast<long long>(n_threads);
    for (unsigned int tid = 0; tid < n_threads; ++tid)
    {
        long long lo = static_cast<long long>(tid) * chunk;
        long long hi = std::min(lo + chunk, n_traj);
        if (lo >= hi)
            break;
        pool.emplace_back(worker, tid, lo, hi);
    }
    for (auto &th : pool)
        th.join();

    // Combine threads and form the maximum-likelihood Poisson extinction time
    // tau = (total time at risk) / (number of extinctions) for each barrier.
    std::vector<long long> died(N_LEVELS, 0);
    std::vector<double> expo(N_LEVELS, 0.0);
    for (unsigned int tid = 0; tid < n_threads; ++tid)
        for (int k = 0; k < N_LEVELS; ++k)
        {
            died[k] += n_died[tid][k];
            expo[k] += exposure[tid][k];
        }

    std::filesystem::path exe_dir = std::filesystem::path(argv[0]).parent_path();
    std::ostringstream adir;
    adir << "alpha_" << std::fixed << std::setprecision(5) << alpha;
    std::filesystem::path out_dir = exe_dir / "outputs" / "extinctionScaling" / adir.str();
    std::filesystem::create_directories(out_dir);
    std::filesystem::path out_path = out_dir / "scaling.tsv";

    std::ofstream file(out_path);
    if (!file.is_open())
    {
        std::cerr << "Failed to open output file: " << out_path.string() << "\n";
        return 1;
    }

    file << "# alpha\t" << std::setprecision(16) << alpha << '\n';
    file << "# t_max\t" << t_max << '\n';
    file << "# n_traj\t" << n_traj << '\n';
    file << "# y0\t" << y0 << '\n';
    file << "# dt\t" << dt << '\n';
    file << "y_ext\tDelta\tn_died\tn_total\ttau\n";
    for (int k = 0; k < N_LEVELS; ++k)
    {
        double Delta = std::log(y0 / y_ext[k]);
        double tau = (died[k] > 0) ? expo[k] / static_cast<double>(died[k])
                                   : std::numeric_limits<double>::quiet_NaN();
        file << std::setprecision(10)
             << y_ext[k] << '\t' << Delta << '\t' << died[k] << '\t' << n_traj << '\t' << tau << '\n';
    }
    file.close();

    std::cout << "Saved " << out_path.string() << " | alpha=" << alpha << "\n";
    return 0;
}
