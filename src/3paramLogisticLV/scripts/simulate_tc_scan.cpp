#include <cmath>
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

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
constexpr int DEFAULT_STEPS_PER_PERIOD = 600;
constexpr int DEFAULT_TRANSIENT_PERIODS = 12;
constexpr int DEFAULT_MEASURE_PERIODS = 6;
constexpr double LOG_Y_MIN = -700.0;
constexpr double LOG_Y_MAX = 700.0;
constexpr double X_MIN = 1e-14;
constexpr double Y_MIN = 1e-14;

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

inline const Env &active_env(double t, double period)
{
    double phase = std::fmod(t, period);
    if (phase < 0.0)
        phase += period;
    return (phase < 0.5 * period) ? ENV_A : ENV_B;
}

inline void rhs(double t, double period, const State &state, double &dx, double &dlogy)
{
    const Env &env = active_env(t, period);
    double y = safe_exp(state.logy);
    dx = env.r * state.x * (1.0 - state.x / env.K) - state.x * y;
    dlogy = state.x - env.gamma;
}

inline void rk4_step(double t, double period, State &state, double dt)
{
    double k1x = 0.0;
    double k1logy = 0.0;
    rhs(t, period, state, k1x, k1logy);

    State s2{state.x + 0.5 * dt * k1x, state.logy + 0.5 * dt * k1logy};
    double k2x = 0.0;
    double k2logy = 0.0;
    rhs(t + 0.5 * dt, period, s2, k2x, k2logy);

    State s3{state.x + 0.5 * dt * k2x, state.logy + 0.5 * dt * k2logy};
    double k3x = 0.0;
    double k3logy = 0.0;
    rhs(t + 0.5 * dt, period, s3, k3x, k3logy);

    State s4{state.x + dt * k3x, state.logy + dt * k3logy};
    double k4x = 0.0;
    double k4logy = 0.0;
    rhs(t + dt, period, s4, k4x, k4logy);

    state.x += (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    state.logy += (dt / 6.0) * (k1logy + 2.0 * k2logy + 2.0 * k3logy + k4logy);

    state.x = clamp_positive(state.x, X_MIN);
    state.logy = std::max(state.logy, std::log(Y_MIN));
}

struct SimulationResult
{
    double growth_rate;
    double x_start_measure;
    double logy_start_measure;
    double x_end;
    double logy_end;
    bool valid;
};

SimulationResult simulate(double x0, double y0, double period, int steps_per_period, int transient_periods, int measure_periods)
{
    const int total_periods = transient_periods + measure_periods;
    const int total_steps = total_periods * steps_per_period;
    const int transient_steps = transient_periods * steps_per_period;
    const int measure_steps = measure_periods * steps_per_period;
    const double dt = period / static_cast<double>(steps_per_period);
    const double measure_window = measure_periods * period;

    State state{x0, std::log(clamp_positive(y0, Y_MIN))};

    double x_start_measure = state.x;
    double logy_start_measure = state.logy;
    bool have_measure_start = false;

    double t = 0.0;
    for (int step = 0; step < total_steps; ++step)
    {
        if (!have_measure_start && step == transient_steps)
        {
            x_start_measure = state.x;
            logy_start_measure = state.logy;
            have_measure_start = true;
        }

        rk4_step(t, period, state, dt);
        t += dt;
    }

    if (!have_measure_start)
    {
        x_start_measure = state.x;
        logy_start_measure = state.logy;
    }

    double growth_rate = (state.logy - logy_start_measure) / measure_window;
    bool valid = std::isfinite(growth_rate) && std::isfinite(state.x) && std::isfinite(state.logy);

    return SimulationResult{growth_rate, x_start_measure, logy_start_measure, state.x, state.logy, valid};
}

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " x0 y0 T task_id [steps_per_period] [transient_periods] [measure_periods]\n";
        return 1;
    }

    const double x0 = std::stod(argv[1]);
    const double y0 = std::stod(argv[2]);
    const double period = std::stod(argv[3]);
    const int task_id = std::stoi(argv[4]);
    const int steps_per_period = (argc > 5) ? std::stoi(argv[5]) : DEFAULT_STEPS_PER_PERIOD;
    const int transient_periods = (argc > 6) ? std::stoi(argv[6]) : DEFAULT_TRANSIENT_PERIODS;
    const int measure_periods = (argc > 7) ? std::stoi(argv[7]) : DEFAULT_MEASURE_PERIODS;

    if (period <= 0.0 || x0 <= 0.0 || y0 <= 0.0 || steps_per_period <= 0 || transient_periods < 0 || measure_periods <= 0)
    {
        std::cerr << "Invalid input parameters. Require x0>0, y0>0, T>0, steps_per_period>0, measure_periods>0.\n";
        return 1;
    }

    std::filesystem::path exe_dir = std::filesystem::path(argv[0]).parent_path();
    std::filesystem::path root_dir = exe_dir;
    std::filesystem::path out_dir = root_dir / "outputs" / "tc_scan" / "raw";
    std::filesystem::create_directories(out_dir);

    std::ostringstream name;
    name << "task_" << std::setfill('0') << std::setw(6) << task_id << ".tsv";
    std::filesystem::path out_path = out_dir / name.str();

    SimulationResult result = simulate(x0, y0, period, steps_per_period, transient_periods, measure_periods);

    std::ofstream file(out_path);
    if (!file.is_open())
    {
        std::cerr << "Failed to open output file: " << out_path.string() << "\n";
        return 1;
    }

    file << "task_id\tx0\ty0\tT\tsteps_per_period\ttransient_periods\tmeasure_periods\tgrowth_rate\tx_start_measure\tlogy_start_measure\tx_end\tlogy_end\tvalid\n";
    file << task_id << '\t'
         << std::setprecision(16) << x0 << '\t'
         << y0 << '\t'
         << period << '\t'
         << steps_per_period << '\t'
         << transient_periods << '\t'
         << measure_periods << '\t'
         << result.growth_rate << '\t'
         << result.x_start_measure << '\t'
         << result.logy_start_measure << '\t'
         << result.x_end << '\t'
         << result.logy_end << '\t'
         << (result.valid ? 1 : 0) << '\n';

    file.close();

    std::cout << "Saved " << out_path.string() << " | growth_rate=" << std::setprecision(10) << result.growth_rate << "\n";
    return 0;
}
