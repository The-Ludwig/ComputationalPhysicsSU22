#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>
#include <yaml-cpp/yaml.h>
#include "NumpySaver.hpp"

constexpr double alpha = 1 / 137.035399; // unitless;
constexpr double hbarc = 197.3269631;    // MeV*fm;
constexpr double c = 299792458;          // m/s

typedef std::complex<double> mattype;
constexpr std::complex<double> I(0, 1);

YAML::Node
get_config(int argc, char const *argv[])
{
    YAML::Node config;
    switch (argc)
    {
    case 2:
        config = YAML::LoadFile(argv[1]);
        break;
    default:
        config = YAML::LoadFile("config.yaml");
        break;
    }
    return config;
}

/**
 * @brief Get the potential at each point in r.
 *
 * This is not yet a good implementation, since I am still learning eigen.
 * There is a lot of optimizing potential.
 *
 * @param r in fm
 * @param R0 in fm
 * @param V0 in MeV
 * @param Z1 number of protons
 * @param Z2 number of protons
 * @return potential in MeV
 */
Eigen::ArrayXd get_potential(Eigen::ArrayXd &r, double V0, double Z1, double Z2)
{
    Eigen::ArrayXd pot(r.rows() + 1);

    pot[0] = alpha * hbarc * Z1 * Z2 / r[0] + V0;

    for (unsigned int i = 1; i < pot.rows() - 1; i++)
    {
        pot[i] = alpha * hbarc * Z1 * Z2 * 2 / (r[i - 1] + r[i]);
    }
    pot[pot.rows() - 1] = 0;

    // for (unsigned int i = 1; i < pot.rows(); i++)
    // {
    //     pot[i] = alpha * hbarc * Z1 * Z2 / r[i - 1];
    // }

    return pot;
}

Eigen::SparseMatrix<mattype> get_matrix(Eigen::ArrayXcd &k, Eigen::ArrayXd &r)
{
    int N = r.size() - 1;
    auto matrix = Eigen::SparseMatrix<mattype>(2 * N + 3, 2 * N + 3);

    // 2N+2 equations, 4 coefficients in each row, +1 A constraint -2 last boundary
    matrix.reserve(4 * (2 * N + 2) + 1 - 2);

    matrix.insert(0, 0) = 1;

    for (int i = 0; i < N; i++)
    {
        auto j = 2 * i + 1;
        matrix.insert(j, j - 1) = std::exp(I * k[i] * r[i]);
        matrix.insert(j, j) = std::exp(-I * k[i] * r[i]);
        matrix.insert(j, j + 1) = -std::exp(I * k[i + 1] * r[i]);
        matrix.insert(j, j + 2) = -std::exp(-I * k[i + 1] * r[i]);

        matrix.insert(j + 1, j - 1) = I * k[i] * std::exp(I * k[i] * r[i]);
        matrix.insert(j + 1, j) = -I * k[i] * std::exp(-I * k[i] * r[i]);
        matrix.insert(j + 1, j + 1) = -I * k[i + 1] * std::exp(I * k[i + 1] * r[i]);
        matrix.insert(j + 1, j + 2) = I * k[i + 1] * std::exp(-I * k[i + 1] * r[i]);
    }

    matrix.insert(2 * N + 1, 2 * N) = std::exp(I * k[N] * r[N]);
    matrix.insert(2 * N + 1, 2 * N + 1) = std::exp(-I * k[N] * r[N]);
    matrix.insert(2 * N + 1, 2 * N + 2) = -std::exp(I * k[N + 1] * r[N]);

    matrix.insert(2 * N + 2, 2 * N) = I * k[N] * std::exp(I * k[N] * r[N]);
    matrix.insert(2 * N + 2, 2 * N + 1) = -I * k[N] * std::exp(-I * k[N] * r[N]);
    matrix.insert(2 * N + 2, 2 * N + 2) = -I * k[N + 1] * std::exp(I * k[N + 1] * r[N]);

    return matrix;
}

Eigen::VectorXcd get_rhs(int N)
{
    Eigen::VectorXcd vec = Eigen::VectorXcd::Zero(2 * N + 3);
    vec[0] = 1;
    return vec;
}

void run_simulation(YAML::Node &simulation_settings, YAML::Node &data)
{
    auto N = simulation_settings["N"].as<int>();
    auto R0 = simulation_settings["R0"].as<double>();
    auto R_last = simulation_settings["R_last"].as<double>();
    auto V0 = simulation_settings["V0"].as<double>();
    auto E_bind_alpha = simulation_settings["E_bind_alpha"].as<double>();
    auto m_neutron = simulation_settings["m_neutron"].as<double>();
    auto m_proton = simulation_settings["m_proton"].as<double>();

    auto name = simulation_settings["Name"].as<std::string>();

    std::cout << "Simulation Setting: " << name << "\n\n";

    for (std::size_t i = 0; i < data.size(); i++)
    {
        auto A_parent = data[i]["A_parent"].as<double>();
        auto Z_parent = data[i]["Z_parent"].as<double>();
        auto E_bind_p = data[i]["E_bind_p"].as<double>();
        auto E_bind_d = data[i]["E_bind_d"].as<double>();

        auto symbol = data[i]["Symbol"].as<std::string>();

        auto E_alpha = E_bind_p - E_bind_d - E_bind_alpha;

        double A_daughter = A_parent - 4;
        double Z_daughter = Z_parent - 2;
        double m_alpha = 3727.379;

        auto R0_this = R0 * std::pow(A_daughter, 1. / 3.);

        // lets first construct the r array as it is constant for every nucleus
        Eigen::ArrayXd r = Eigen::ArrayXd::LinSpaced(N + 1, R0_this, R_last);

        Eigen::ArrayXcd pot = get_potential(r, V0, Z_daughter, 2); // MeV
        Eigen::ArrayXcd k = (2 * m_alpha * (E_alpha - pot)).sqrt() / hbarc;

        NumpySaver(std::string("build/output/" + name + "-" + symbol + "-r.npy")) << r;
        NumpySaver(std::string("build/output/" + name + "-" + symbol + "-pot.npy")) << pot.real();

        auto mat = get_matrix(k, r);
        auto b = get_rhs(N);

        Eigen::SparseLU<Eigen::SparseMatrix<mattype>, Eigen::COLAMDOrdering<int>> solver;
        // Compute the ordering permutation vector from the structural pattern of A
        solver.analyzePattern(mat);
        // Compute the numerical factorization
        solver.factorize(mat);
        // Use the factors to solve the linear system
        Eigen::VectorXcd x = solver.solve(b);

        // R and T factor
        double R = std::norm(x(1));
        double T = std::norm(x(Eigen::last)) * std::abs(k(Eigen::last)) / std::abs(k[0]);
        // check if R+T=1
        if (std::abs(R + T - 1) > 1e-3)
        {
            std::cerr << "Warning! R+T not equal to 1; R+T = " << R + T << "\n";
        }

        double v = std::sqrt(2 * E_alpha / m_alpha) * c;              // m/s
        double half_life = std::log(2) * 2 * R0_this / T / v * 1e-15; // s

        std::cerr << "Lifetime of " << symbol << " " << half_life << "s" << std::endl;
    }
}

int main(int argc, char const *argv[])
{
    auto config = get_config(argc, argv);
    auto data = config["data"];

    for (std::size_t i = 0; i < config["simulation_settings"].size(); i++)
    {
        std::cout << "##################################\n"
                  << "# Running Simulation " << i + 1 << "/" << config["simulation_settings"].size()
                  << "\n##################################\n";
        auto settings = config["simulation_settings"][i];
        run_simulation(settings, data);
        std::cout << "\n\n";
    }

    return 0;
}
