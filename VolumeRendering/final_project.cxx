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

template<typename T>
void UnitNorm(const T *a, const T *b, T *c)
{
    vtkMath::Cross(a,b,c);
    vtkMath::Normalize(c);
}

struct Camera
{
    double          near, far;
    double          angle;
    double          W,H;
    double          position[3];
    double          focus[3];
    double          up[3];
    double          look[3];
    double          u[3];
    double          v[3];
    double          fov[2];
    double          dx[3];
    double          dy[3];
    void            GetDelta(double*, double*, double*, double*);
    void            Pixel2Ray();
};

void
Camera::Pixel2Ray()
{
    for( size_t i = 0; i < 3; ++i)
        look[i] = this->focus[i] - this->position[i];
    UnitNorm(this->look, this->up, this->u);
    UnitNorm(this->look, this->u, this->v);
    for( size_t i = 0; i < 3; ++i )
    {
        this->dx[i] = 2.*tan(fov[0]/2)/this->W * u[i];
        this->dy[i] = 2.*tan(fov[1]/2)/this->H * v[i];
    }
    //printf("dx = <%f,%f,%f>\n",dx[0],dx[1],dx[2]);
    //printf("dy = <%f,%f,%f>\n",dy[0],dy[1],dy[2]);
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
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;
    // my additions
    rv.fov[0] = rv.angle * PI/180; // Radians
    rv.fov[1] = rv.angle * PI/180; // Radians
    rv.W = 1024;
    rv.H = 1024;

    return rv;
}

int main()
{
    vtkSmartPointer<vtkDataSetReader> rdr =
        vtkSmartPointer<vtkDataSetReader>::New();
    rdr->SetFileName("astro64.vtk");
    rdr->Update();

    int dims[3];
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
    camera.Pixel2Ray();



    return 0;
}
