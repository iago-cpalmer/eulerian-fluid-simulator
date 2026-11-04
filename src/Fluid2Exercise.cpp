#include "Scene.h"

#include "Numeric/PCGSolver.h"

namespace asa
{
namespace
{

template <class T>
T bilerp(const Array2<T> &grid, float x, float y)
{
    int x0 = (int)floor(x);
    int y0 = (int)floor(y);

    int x1 = x0 + 1;
    int y1 = y0 + 1;

    float sx = x - x0;
    float sy = y - y0;

    x0 = std::max(0, std::min(x0, (int)grid.getSize().x - 1));
    x1 = std::max(0, std::min(x1, (int)grid.getSize().x - 1));
    y0 = std::max(0, std::min(y0, (int)grid.getSize().y - 1));
    y1 = std::max(0, std::min(y1, (int)grid.getSize().y - 1));

    const T &v00 = grid.getValue(x0, y0);
    const T &v10 = grid.getValue(x1, y0);
    const T &v01 = grid.getValue(x0, y1);
    const T &v11 = grid.getValue(x1, y1);

    return v00 * ((1 - sx) * (1 - sy)) + v10 * (sx * (1 - sy)) + v01 * ((1 - sx) * sy) + v11 * (sx * sy);
}

// Velocity at a world-space position. Each component lives on its own staggered
// faces, so each one is looked up with its own world-to-index conversion.
Vector2 sampleVelocity(const Fluid2 &f, const Vector2 &pos)
{
    const Vector2 idU = f.grid.getFaceIndexX(pos);
    const Vector2 idV = f.grid.getFaceIndexY(pos);
    return Vector2(bilerp(f.velocityX, idU.x, idU.y), bilerp(f.velocityY, idV.x, idV.y));
}


void applySolidWalls(Array2<float> &u, Array2<float> &v)
{
    int nx = u.getSize().x;
    int ny = v.getSize().y;

    // LEFT / RIGHT (u)
    for (int j = 0; j < u.getSize().y; j++) {
        u.setValue(0, j, 0.0f);
        u.setValue(nx - 1, j, 0.0f);
    }

    // BOTTOM / TOP (v)
    for (int i = 0; i < v.getSize().x; i++) {
        v.setValue(i, 0, 0.0f);
        v.setValue(i, ny - 1, 0.0f);
    }
}



}  // namespace

////////////////////////////////////////////////////////////
////////////////////// ADVECTION ///////////////////////////
////////////////////////////////////////////////////////////

void Fluid2::fluidAdvection(const float dt)
{
    {
        Array2<Vector3> newInk(inkRGB.getSize());

        for (uint j = 0; j < inkRGB.getSize().y; j++)
            for (uint i = 0; i < inkRGB.getSize().x; i++) {
                const Vector2 pos = grid.getCellPos(Index2(i, j));
                const Vector2 prev = pos - dt * sampleVelocity(*this, pos);
                const Vector2 id = grid.getCellIndex(prev);
                newInk.setValue(i, j, bilerp(inkRGB, id.x, id.y));
            }

        inkRGB = newInk;
    }

    {
        Array2<float> newVX(velocityX.getSize());
        Array2<float> newVY(velocityY.getSize());

        for (uint j = 0; j < velocityX.getSize().y; j++)
            for (uint i = 0; i < velocityX.getSize().x; i++) {
                const Vector2 pos = grid.getFacePosX(Index2(i, j));
                const Vector2 prev = pos - dt * sampleVelocity(*this, pos);
                const Vector2 id = grid.getFaceIndexX(prev);
                newVX.setValue(i, j, bilerp(velocityX, id.x, id.y));
            }

        for (uint j = 0; j < velocityY.getSize().y; j++)
            for (uint i = 0; i < velocityY.getSize().x; i++) {
                const Vector2 pos = grid.getFacePosY(Index2(i, j));
                const Vector2 prev = pos - dt * sampleVelocity(*this, pos);
                const Vector2 id = grid.getFaceIndexY(prev);
                newVY.setValue(i, j, bilerp(velocityY, id.x, id.y));
            }

        velocityX = newVX;
        velocityY = newVY;
    }
}

////////////////////////////////////////////////////////////
////////////////////// EMISSION ////////////////////////////
////////////////////////////////////////////////////////////

void Fluid2::fluidEmission()
{
    if (Scene::testcase >= Scene::SMOKE) {
        Index2 size = grid.getSize();

        Vector2 domainMin = grid.getDomain().minPosition;
        Vector2 domainMax = grid.getDomain().maxPosition;

        auto toGrid = [&](float x, float y, int &i, int &j) {
            i = (int)((x - domainMin.x) / (domainMax.x - domainMin.x) * size.x);
            j = (int)((y - domainMin.y) / (domainMax.y - domainMin.y) * size.y);
        };

        int i0, j0, i1, j1;

        // Emisor 1
        toGrid(-0.1f, -1.9f, i0, j0);
        toGrid(0.1f, -1.75f, i1, j1);

        for (int j = j0; j <= j1; j++)
            for (int i = i0; i <= i1; i++)
                if (i >= 0 && i < size.x && j >= 0 && j < size.y) {
                    inkRGB.setValue(i, j, Vector3(1.0f, 1.0f, 0.0f));
                    velocityX.setValue(i, j, 0.0f);
                    velocityY.setValue(i, j, 8.0f);
                }

        // Emisor 2
        toGrid(-0.22f, -1.9f, i0, j0);
        toGrid(-0.1f, -1.75f, i1, j1);

        for (int j = j0; j <= j1; j++)
            for (int i = i0; i <= i1; i++)
                if (i >= 0 && i < size.x && j >= 0 && j < size.y) {
                    inkRGB.setValue(i, j, Vector3(1.0f, 0.0f, 1.0f));
                    velocityX.setValue(i, j, 0.0f);
                    velocityY.setValue(i, j, 8.0f);
                }

        // Emisor 3
        toGrid(0.1f, -1.9f, i0, j0);
        toGrid(0.2f, -1.75f, i1, j1);

        for (int j = j0; j <= j1; j++)
            for (int i = i0; i <= i1; i++)
                if (i >= 0 && i < size.x && j >= 0 && j < size.y) {
                    inkRGB.setValue(i, j, Vector3(0.0f, 1.0f, 1.0f));
                    velocityX.setValue(i, j, 0.0f);
                    velocityY.setValue(i, j, 8.0f);
                }
    }
}

////////////////////////////////////////////////////////////
////////////////////// FORCES //////////////////////////////
////////////////////////////////////////////////////////////

void Fluid2::fluidVolumeForces(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        for (uint j = 0; j < velocityY.getSize().y; j++)
            for (uint i = 0; i < velocityY.getSize().x; i++)
                velocityY.setValue(i, j, velocityY.getValue(i, j) + Scene::kGravity * dt);
    }
}

////////////////////////////////////////////////////////////
////////////////////// VISCOSITY ///////////////////////////
////////////////////////////////////////////////////////////

void Fluid2::fluidViscosity(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        // Start from a copy: the loops below skip the outermost faces, which must keep their values.
        Array2<float> newVX(velocityX);
        Array2<float> newVY(velocityY);

        const Vector2 dx = grid.getDx();
        const float invDx2 = 1.0f / (dx.x * dx.x);
        const float invDy2 = 1.0f / (dx.y * dx.y);

        for (uint j = 1; j < velocityX.getSize().y - 1; j++)
            for (uint i = 1; i < velocityX.getSize().x - 1; i++) {
                const float c = velocityX.getValue(i, j);
                const float lap = (velocityX.getValue(i - 1, j) - 2 * c + velocityX.getValue(i + 1, j)) * invDx2 +
                                  (velocityX.getValue(i, j - 1) - 2 * c + velocityX.getValue(i, j + 1)) * invDy2;
                newVX.setValue(i, j, c + Scene::kViscosity * lap * dt / Scene::kDensity);
            }

        for (uint j = 1; j < velocityY.getSize().y - 1; j++)
            for (uint i = 1; i < velocityY.getSize().x - 1; i++) {
                const float c = velocityY.getValue(i, j);
                const float lap = (velocityY.getValue(i - 1, j) - 2 * c + velocityY.getValue(i + 1, j)) * invDx2 +
                                  (velocityY.getValue(i, j - 1) - 2 * c + velocityY.getValue(i, j + 1)) * invDy2;
                newVY.setValue(i, j, c + Scene::kViscosity * lap * dt / Scene::kDensity);
            }

        velocityX = newVX;
        velocityY = newVY;
    }
}

////////////////////////////////////////////////////////////
//////////////////// PRESSURE PROJECTION ///////////////////
////////////////////////////////////////////////////////////

void Fluid2::fluidPressureProjection(float dt)
{
    if (Scene::testcase < Scene::SMOKE)
        return;

    int nx = grid.getSize().x;
    int ny = grid.getSize().y;
    int N = nx * ny;

    applySolidWalls(velocityX, velocityY);

    std::vector<float> rhs(N, 0.0f);

    auto idx = [&](int i, int j) { return j * nx + i; };

    const Vector2 dx = grid.getDx();
    const float invDx2 = 1.0f / (dx.x * dx.x);
    const float invDy2 = 1.0f / (dx.y * dx.y);

    // divergence
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            float div = (velocityX.getValue(i + 1, j) - velocityX.getValue(i, j)) / dx.x +
                        (velocityY.getValue(i, j + 1) - velocityY.getValue(i, j)) / dx.y;

            rhs[idx(i, j)] = -Scene::kDensity / dt * div;
        }

    SparseMatrix<float> A(N, 5);

    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            int id = idx(i, j);

            float diag = 0.0f;

            if (i > 0) {
                A.add_to_element(id, idx(i - 1, j), -invDx2);
                diag += invDx2;
            }
            if (i < nx - 1) {
                A.add_to_element(id, idx(i + 1, j), -invDx2);
                diag += invDx2;
            }
            if (j > 0) {
                A.add_to_element(id, idx(i, j - 1), -invDy2);
                diag += invDy2;
            }
            if (j < ny - 1) {
                A.add_to_element(id, idx(i, j + 1), -invDy2);
                diag += invDy2;
            }

            A.add_to_element(id, id, diag);
        }

    // solve
    PCGSolver<float> solver;
    solver.set_solver_parameters(1e-5f, 200, 0.97f, 0.25f);

    std::vector<float> p(N, 0.0f);
    float res;
    int iters;

    solver.solve(A, rhs, p, res, iters);

    // store pressure
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            pressure.setValue(i, j, p[idx(i, j)]);

    float scale = dt / Scene::kDensity;

    // subtract gradient
    for (int j = 0; j < ny; j++)
        for (int i = 1; i < nx; i++) {
            velocityX.setValue(
                i, j, velocityX.getValue(i, j) - scale * (pressure.getValue(i, j) - pressure.getValue(i - 1, j)) / dx.x);
        }

    for (int j = 1; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            velocityY.setValue(
                i, j, velocityY.getValue(i, j) - scale * (pressure.getValue(i, j) - pressure.getValue(i, j - 1)) / dx.y);
        }

    applySolidWalls(velocityX, velocityY);
}
/*
void Fluid2::fluidPressureProjection(const float dt)
{

    if (Scene::testcase >= Scene::SMOKE) {
        Index2 size = grid.getSize();
        Array2<float> divergence(size);

        applySolidWalls(velocityX, velocityY);


        for (uint j = 0; j < size.y; j++)
            for (uint i = 0; i < size.x; i++) {
                float div = velocityX.getValue(i + 1, j) - velocityX.getValue(i, j) + velocityY.getValue(i, j + 1) -
                            velocityY.getValue(i, j);
                divergence.setValue(i, j, div * Scene::kDensity);
            }

        SparseMatrix<float> A(size.x * size.y, 5);
        for (int iter = 0; iter < 40; iter++) {
            Array2<float> newPressure = pressure;
            for (uint j = 1; j < size.y - 1; j++)
                for (uint i = 1; i < size.x - 1; i++) {
                    float p = pressure.getValue(i - 1, j) + pressure.getValue(i + 1, j) + pressure.getValue(i, j - 1) +
                              pressure.getValue(i, j + 1) - divergence.getValue(i, j);
                    newPressure.setValue(i, j, p * 0.25f);
                }
            pressure = newPressure;
        }
        for (uint j = 0; j < size.y; j++)
            for (uint i = 1; i < velocityX.getSize().x - 1; i++) {
                float grad = pressure.getValue(i, j) - pressure.getValue(i - 1, j);
                velocityX.setValue(i, j, velocityX.getValue(i, j) - grad);
            }
        for (uint j = 1; j < velocityY.getSize().y - 1; j++)
            for (uint i = 0; i < size.x; i++) {
                float grad = pressure.getValue(i, j) - pressure.getValue(i, j - 1);
                velocityY.setValue(i, j, velocityY.getValue(i, j) - grad);
            }
    }
}*/

}  // namespace asa