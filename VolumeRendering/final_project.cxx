#include <iostream>
// VTK
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
// My VTK
#include <vtkMath.h>
#include <vtkDataSetReader.h>
#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkType.h>

#define PI 3.141592654
#define DIM 3

//template<typename T>
void CrossUnitNorm(const double *a, const double *b, double *c)
{
    vtkMath::Cross(a,b,c);
    vtkMath::Normalize(c);
}

class Ray
{
    private:
        double          direction[DIM];
        double          origin[DIM];
        double          dx[DIM];
        double          dy[DIM];
        int             samples;

    public:
        Ray(double *origin, int samples);
        void            Setdx(double*);
        void            Setdy(double*);
        void            SetDirection(double*);
        inline double   GetRayDirection(int i){return direction[i];};
        void            PrintRayInfo();
};

Ray::Ray(double *o, int s)
{
    for (std::size_t i = 0; i < DIM; ++i)
        origin[i] = o[i];
    samples = s;
}

void
Ray::Setdx(double *x)
{
    for (std::size_t i = 0; i < DIM; ++i)
        dx[i] = x[i];
}

void
Ray::Setdy(double *y)
{
    for (std::size_t i = 0; i < DIM; ++i)
        dy[i] = y[i];
}

void
Ray::SetDirection(double *d)
{
    for (std::size_t i = 0; i < DIM; ++i)
        direction[i] = d[i];
}

void
Ray::PrintRayInfo()
{
    printf("========== Ray Info ==========\n");
    printf("Direction <%f,%f,%f>\n", direction[0], direction[1], direction[2]);
    printf("Origin <%f,%f,%f>\n", origin[0], origin[1], origin[2]);
    printf("delta_x <%f,%f,%f>\n", dx[0], dx[1], dx[2]);
    printf("delta_y <%f,%f,%f>\n", dy[0], dy[1], dy[2]);
    printf("samples %d\n", samples);
    printf("==============================\n");
}


struct Camera
{
    double          near, far;
    double          angle; // fov_x and fov_y
    double          W,H;
    double          position[DIM];
    double          focus[DIM];
    double          up[DIM];
    double          look[DIM];
    double          u[DIM];
    double          v[DIM];
    void            GetDelta(double*, double*, double*, double*);
    void            Pixel2Ray(Ray&, int, int);
    inline double*  GetOrigin(){return position;};
};

void
Camera::Pixel2Ray(Ray &ray, int i, int j)
{
    double dx[DIM];
    double dy[DIM];
    double direction[DIM];
    for( size_t i = 0; i < DIM; ++i)
        look[i] = focus[i] - position[i];
    CrossUnitNorm(look, up, u);
    //printf("ru = <%f,%f,%f>\n", u[0], u[1], u[2]);
    CrossUnitNorm(look, u, v);
    //printf("rv = <%f,%f,%f>\n", v[0], v[1], v[2]);
    // NOTE: look is normalized at this point
    for( size_t i = 0; i < DIM; ++i )
    {
        dx[i] = 2.*tan(angle/2)/W * u[i];
        dy[i] = 2.*tan(angle/2)/H * v[i];
    }
    ray.Setdx(dx);
    ray.Setdy(dy);
    //printf("dx = <%f,%f,%f>\n",dx[0],dx[1],dx[2]);
    //printf("dy = <%f,%f,%f>\n",dy[0],dy[1],dy[2]);
    // Direction
    double look_norm = vtkMath::Normalize(look);
    //printf("look = <%f,%f,%f>\n", look[0], look[1], look[2]);
    //printf("look norm = %f\n", vtkMath::Normalize(look));
    direction[0] = look[0] + ((2*i+1 - W)/2)*dx[0] + ((2*j+1-H)/2)*dy[0];
    direction[1] = look[1] + ((2*i+1 - W)/2)*dx[1] + ((2*j+1-H)/2)*dy[1];
    direction[2] = look[2] + ((2*i+1 - W)/2)*dx[2] + ((2*j+1-H)/2)*dy[2];
    ray.SetDirection(direction);
    //printf("look = <%f,%f,%f>, W = %f, H = %f\n", look[0], look[1], look[2], W, H);
    //printf("Direction <%f,%f,%f>\n", direction[0], direction[1], direction[2]);
}

struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    double         *opacities; // size is numBins

    // Take in a value and applies the transfer function.
    // Step #1: figure out which bin "value" lies in.
    // If "min" is 2 and "max" is 4, and there are 10 bins, then
    //   bin 0 = 2->2.2
    //   bin 1 = 2.2->2.4
    //   bin 2 = 2.4->2.6
    //   bin 3 = 2.6->2.8
    //   bin 4 = 2.8->3.0
    //   bin 5 = 3.0->3.2
    //   bin 6 = 3.2->3.4
    //   bin 7 = 3.4->3.6
    //   bin 8 = 3.6->3.8
    //   bin 9 = 3.8->4.0
    // and, for example, a "value" of 3.15 would return the color in bin 5
    // and the opacity at "opacities[5]".
    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
    {
        //int bin = GetBin(value);
        //RGB[0] = colors[3*bin+0];
        //RGB[1] = colors[3*bin+1];
        //RGB[2] = colors[3*bin+2];
        //opacity = opacities[bin];
    }
};

Camera
SetupCamera(void)
{
    Camera rv;
    rv.focus[0] = 0;
    rv.focus[1] = 0;
    rv.focus[2] = 0;
    rv.up[0] = 0;
    rv.up[1] = -1;
    rv.up[2] = 0;
    rv.angle = 30 * PI/180; // Fix to be radians
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;
    // my additions
    //rv.fov[0] = rv.angle * PI/180; // Radians
    //rv.fov[1] = rv.angle * PI/180; // Radians
    //rv.W = 1024;
    //rv.H = 1024;
    rv.W = 100;
    rv.H = 100;

    return rv;
}

int main()
{
    int samples = 256;
    vtkSmartPointer<vtkDataSetReader> rdr =
        vtkSmartPointer<vtkDataSetReader>::New();
    rdr->SetFileName("astro64.vtk");
    rdr->Update();

    int dims[DIM];
    //vtkSmartPointer<vtkRectilinearGrid> rgrid = 
    //    vtkSmartPointer<vtkRectilinearGrid>::Take(vtkRectilinearGrid::SafeDownCast(rdr->GetOutput()));
    //rgrid->GetDimensions(dims);
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float*) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float*) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float*) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float*) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    Camera camera = SetupCamera();
    Ray ray(camera.GetOrigin(), samples);// = new Ray();
    camera.Pixel2Ray(ray, 50, 50);
    //printf("ray direction <%f,%f,%f>\n", ray.GetRayDirection(0), ray.GetRayDirection(1), ray.GetRayDirection(2));
    ray.PrintRayInfo();



    return 0;
}
