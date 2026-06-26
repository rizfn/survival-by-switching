#include <atomic>
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

// Standard parameters: the switching paradox is present, with a finite
// critical switching rate alpha_c ~ 1.556 (see math/stochastic.md).
constexpr Env ENV_A{1.0, 1.0, 1.0};
constexpr Env ENV_B{4.0, 2.0, 2.2};

constexpr double DEFAULT_T_MAX = 6000.0;     // horizon over which we ask "did y survive?"
constexpr long long DEFAULT_N_TRAJ = 20000; // ensemble size
constexpr double DEFAULT_DT = 0.01;         // RK4 step inside a dwell
constexpr double DEFAULT_Y_EXT = 1e-3;      // extinction threshold on y
constexpr double LOG_Y_MIN = -700.0;
constexpr double LOG_Y_MAX = 700.0;
constexpr double X_MIN = 1e-14;

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

// Right-hand side for a fixed environment (the environment is constant within a
// single dwell, so it does not need to be re-selected from t here). We integrate
// the predator in log-space, dlog(y)/dt = x - gamma, to keep the huge dynamic
// range of y under control.
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

// Simulate one trajectory with Poisson (rate alpha) environment switching.
// Returns the time at which y first falls below y_ext (extinction), or -1.0 if
// the predator survives the whole horizon t_max.
double simulate_trajectory(double alpha, double x0, double y0, double t_max,
                           double dt, double log_y_ext, std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::exponential_distribution<double> dwell(alpha);
    std::bernoulli_distribution coin(0.5);

    State state{x0, std::log(y0)};
    const Env *env = coin(rng) ? &ENV_A : &ENV_B;

    double t = 0.0;
    while (t < t_max)
    {
        // Length of this dwell, clipped to the remaining horizon.
        double dwell_len = dwell(rng);
        double dwell_end = std::min(t + dwell_len, t_max);

        while (t < dwell_end)
        {
            double step = std::min(dt, dwell_end - t);
            rk4_step(*env, state, step);
            t += step;
            if (state.logy <= log_y_ext)
                return t; // y has died out
        }

        // Flip the environment for the next dwell.
        env = (env == &ENV_A) ? &ENV_B : &ENV_A;
    }

    return -1.0; // survived to t_max
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0]
                  << " alpha sim_id [t_max] [n_traj] [dt] [y_ext] [seed]\n";
        return 1;
    }

    const double alpha = std::stod(argv[1]);
    const int sim_id = std::stoi(argv[2]);
    const double t_max = (argc > 3) ? std::stod(argv[3]) : DEFAULT_T_MAX;
    const long long n_traj = (argc > 4) ? std::stoll(argv[4]) : DEFAULT_N_TRAJ;
    const double dt = (argc > 5) ? std::stod(argv[5]) : DEFAULT_DT;
    const double y_ext = (argc > 6) ? std::stod(argv[6]) : DEFAULT_Y_EXT;
    const std::uint64_t base_seed =
        (argc > 7) ? static_cast<std::uint64_t>(std::stoull(argv[7])) : 0x9E3779B97F4A7C15ULL;

    if (alpha <= 0.0 || t_max <= 0.0 || n_traj <= 0 || dt <= 0.0 || y_ext <= 0.0)
    {
        std::cerr << "Invalid parameters. Require alpha>0, t_max>0, n_traj>0, dt>0, y_ext>0.\n";
        return 1;
    }

    // Start both populations at the midpoint of the two carrying capacities, so
    // the predator begins reasonably far from either equilibrium and from zero.
    const double x0 = 0.5 * (ENV_A.K + ENV_B.K);
    const double y0 = 0.5 * (ENV_A.K + ENV_B.K);
    const double log_y_ext = std::log(y_ext);

    // Ensemble in parallel, leaving 2 cores free for the rest of the machine.
    unsigned int hw = std::thread::hardware_concurrency();
    unsigned int n_threads = (hw > 3) ? hw - 2 : 1;
    if (static_cast<long long>(n_threads) > n_traj)
        n_threads = static_cast<unsigned int>(n_traj);

    std::vector<double> death_time(static_cast<std::size_t>(n_traj), -1.0);

    // Each trajectory has its own deterministic seed, so the result is
    // reproducible regardless of how the work is split across threads.
    auto worker = [&](long long lo, long long hi)
    {
        for (long long i = lo; i < hi; ++i)
        {
            std::uint64_t seed = base_seed + 0x100000001B3ULL * static_cast<std::uint64_t>(i);
            death_time[static_cast<std::size_t>(i)] =
                simulate_trajectory(alpha, x0, y0, t_max, dt, log_y_ext, seed);
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
        pool.emplace_back(worker, lo, hi);
    }
    for (auto &th : pool)
        th.join();

    // Output: compact, one death time per trajectory (-1 = survived to t_max).
    std::filesystem::path exe_dir = std::filesystem::path(argv[0]).parent_path();
    std::ostringstream pdir;
    pdir << "alpha_" << std::fixed << std::setprecision(5) << alpha
         << "_tmax_" << std::defaultfloat << t_max
         << "_ntraj_" << n_traj
         << "_yext_" << std::defaultfloat << y_ext;
    std::filesystem::path out_dir = exe_dir / "outputs" / "survivalProb" / pdir.str();
    std::filesystem::create_directories(out_dir);

    std::ostringstream name;
    name << "sim_" << sim_id << ".tsv";
    std::filesystem::path out_path = out_dir / name.str();

    std::ofstream file(out_path);
    if (!file.is_open())
    {
        std::cerr << "Failed to open output file: " << out_path.string() << "\n";
        return 1;
    }

    file << "# alpha\t" << std::setprecision(16) << alpha << '\n';
    file << "# t_max\t" << t_max << '\n';
    file << "# n_traj\t" << n_traj << '\n';
    file << "# x0\t" << x0 << '\n';
    file << "# y0\t" << y0 << '\n';
    file << "# y_ext\t" << y_ext << '\n';
    file << "# dt\t" << dt << '\n';
    file << "# death_time (-1 = survived to t_max)\n";
    for (long long i = 0; i < n_traj; ++i)
        file << death_time[static_cast<std::size_t>(i)] << '\n';
    file.close();

    long long n_died = 0;
    for (double d : death_time)
        if (d >= 0.0)
            ++n_died;
    std::cout << "Saved " << out_path.string()
              << " | alpha=" << alpha
              << " survived=" << (n_traj - n_died) << "/" << n_traj
              << " (" << std::setprecision(4) << 100.0 * (n_traj - n_died) / n_traj << "%)\n";
    return 0;
}
