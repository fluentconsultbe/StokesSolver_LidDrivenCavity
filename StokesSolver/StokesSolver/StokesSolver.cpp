#include <vector>
#include <cmath>
#include <iostream>

struct Field
{
    int nx, ny;
    double dx, dy;
    std::vector<double> data;

    Field(int nx_, int ny_, double Lx, double Ly)
        : nx(nx_), ny(ny_),
        dx(Lx / (nx_ - 1)),
        dy(Ly / (ny_ - 1)),
        data(nx_* ny_, 0.0)
    {
    }

    double& operator()(int i, int j) { return data[j * nx + i]; }
    const double& operator()(int i, int j) const { return data[j * nx + i]; }
};

void applyBC(Field& u, Field& v)
{
    int nx = u.nx, ny = u.ny;

    for (int i = 0; i < nx; i++)
    {
        u(i, ny - 1) = 1.0;
        v(i, ny - 1) = 0.0;
    }

    for (int i = 0; i < nx; i++)
    {
        u(i, 0) = 0.0;
        v(i, 0) = 0.0;
    }

    for (int j = 0; j < ny; j++)
    {
        u(0, j) = 0.0;
        v(0, j) = 0.0;
        u(nx - 1, j) = 0.0;
        v(nx - 1, j) = 0.0;
    }
}

double laplace(const Field& f, int i, int j)
{
    double h = f.dx;
    return (
        f(i + 1, j) +
        f(i - 1, j) +
        f(i, j + 1) +
        f(i, j - 1) -
        4.0 * f(i, j)
        ) / (h * h);
}

void normalizePressure(Field& p)
{
    double sum = 0.0;
    for (double val : p.data) sum += val;
    double mean = sum / p.data.size();
    for (double& val : p.data) val -= mean;
}

double divergence(const Field& u, const Field& v, int i, int j)
{
    double dx = u.dx;
    double dy = u.dy;

    double du_dx = (u(i + 1, j) - u(i - 1, j)) / (2.0 * dx);
    double dv_dy = (v(i, j + 1) - v(i, j - 1)) / (2.0 * dy);

    return du_dx + dv_dy;
}

double dpdx(const Field& p, int i, int j)
{
    return (p(i + 1, j) - p(i - 1, j)) / (2.0 * p.dx);
}

double dpdy(const Field& p, int i, int j)
{
    return (p(i, j + 1) - p(i, j - 1)) / (2.0 * p.dy);
}

void solvePressurePoisson(
    Field& p,
    const Field& uStar,
    const Field& vStar,
    int iters,
    double dt
)
{
    Field pNew = p;

    for (int it = 0; it < iters; ++it)
    {
        for (int j = 1; j < p.ny - 1; ++j)
        {
            for (int i = 1; i < p.nx - 1; ++i)
            {
                double rhs = (1.0 / dt) * divergence(uStar, vStar, i, j);
                double h = p.dx;

                pNew(i, j) = 0.25 * (
                    p(i + 1, j) + p(i - 1, j) +
                    p(i, j + 1) + p(i, j - 1) -
                    h * h * rhs
                    );
            }
        }

        for (int i = 0; i < p.nx; ++i)
        {
            pNew(i, 0) = pNew(i, 1);
            pNew(i, p.ny - 1) = pNew(i, p.ny - 2);
        }
        for (int j = 0; j < p.ny; ++j)
        {
            pNew(0, j) = pNew(1, j);
            pNew(p.nx - 1, j) = pNew(p.nx - 2, j);
        }

        p.data.swap(pNew.data);
        normalizePressure(p);
    }
}

double computeDivergenceResidual(const Field& u, const Field& v)
{
    double maxDiv = 0.0;
    for (int j = 1; j < u.ny - 1; ++j)
    {
        for (int i = 1; i < u.nx - 1; ++i)
        {
            double div = std::abs(divergence(u, v, i, j));
            if (div > maxDiv) maxDiv = div;
        }
    }
    return maxDiv;
}

int main()
{
    int nx = 64, ny = 64;
    double Lx = 1.0, Ly = 1.0;

    double Re = 100.0;
    double nu = 1.0 / Re;   // 0.01
    double dt = 1e-3;

    Field u(nx, ny, Lx, Ly);
    Field v(nx, ny, Lx, Ly);
    Field p(nx, ny, Lx, Ly);
    Field uStar(nx, ny, Lx, Ly);
    Field vStar(nx, ny, Lx, Ly);

    applyBC(u, v);

    int maxIter = 5000;
    for (int iter = 0; iter < maxIter; ++iter)
    {
        applyBC(u, v);

        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                double lap_u = laplace(u, i, j);
                double lap_v = laplace(v, i, j);

                uStar(i, j) = u(i, j) + dt * nu * lap_u;
                vStar(i, j) = v(i, j) + dt * nu * lap_v;
            }
        }

        applyBC(uStar, vStar);

        solvePressurePoisson(p, uStar, vStar, 200, dt);

        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                u(i, j) = uStar(i, j) - dt * dpdx(p, i, j);
                v(i, j) = vStar(i, j) - dt * dpdy(p, i, j);
            }
        }

        applyBC(u, v);

        if (iter % 100 == 0)
        {
            double R = computeDivergenceResidual(u, v);
            std::cout << "Iter: " << iter
                << "   Max divergence = " << R << std::endl;
        }
    }

    std::cout << "Done.\n";
    return 0;
}
