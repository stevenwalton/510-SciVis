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
//#include <vtkDoubleArray.h>
// My VTK
#include <vtkMath.h>
#include <vtkDataSetReader.h>
#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkType.h>

#define PI 3.141592654
#define DIM 3
#define WIDTH 100
#define HEIGHT 100

/*
int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
}
*/

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
}

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);
}

/*
LinearInterpolation(const float a, const float b, const float x,
                    const float F_a, const float F_b)
{
    float t = (x-a)/(b-a);
    return F_a + t*(F_b - F_a);
}
*/

int
BinarySearch(const float pt, const int lower, const int upper,
             const float *arr)
{
    if (lower > upper)
    {
        printf("Error in binary search! %d > %d", lower, upper);
        exit(1);
    }
    int index = (lower + upper)/2;
    if (lower == upper)
        return index;
    else if (pt > arr[index])
    {
        if (pt < arr[index+1])
        {
            return index;
        }
        index = BinarySearch(pt, index+1, upper, arr);
    }
    else if (pt < arr[index])
    {
        index = BinarySearch(pt, lower, index-1, arr);
    }
    return index;
}

//template<typename T>
void CrossUnitNorm(const float *a, const float *b, float *c)
{
    vtkMath::Cross(a,b,c);
    vtkMath::Normalize(c);
}


struct TransferFunction {
    float          min;
    float          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    float         *opacities; // size is numBins
    float          scale;

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
    void ApplyTransferFunction(float value, unsigned char *RGB, double &opacity)
    {
        if (value < min || value > max) return;
        //int bin = GetBin(value);
        int bin = (value - min);
        RGB[0] = colors[3*bin+0];
        RGB[1] = colors[3*bin+1];
        RGB[2] = colors[3*bin+2];
        opacity = opacities[bin];
    }
    /*
    int GetBin(float val)
    {
        if (val < min) return min;
        else if (val > max) return max;
        else return (val - min)*scale;
    }
    */
};

TransferFunction
SetupTransferFunction(void)
{
    int  i;

    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    //rv.opacities = new double[256];
    rv.opacities = new float[256];
    unsigned char charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = charOpacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = { 
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76 
       };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0 ; i < numControlPoints-1 ; i++)
    {
        int start = controlPointPositions[i]*rv.numBins;
        int end   = controlPointPositions[i+1]*rv.numBins+1;
        //cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << endl;
        if (end >= rv.numBins)
            end = rv.numBins-1;
        for (int j = start ; j <= end ; j++)
        {
            double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (int k = 0 ; k < 3 ; k++)
                rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                                 + controlPointColors[3*i+k];
        }
    }    

    rv.scale = 256/(rv.max - rv.min); // Check
    return rv;
}

class Ray
{
    private:
        float          dx[DIM];
        float          dy[DIM];

    public:
        Ray(float *origin, int samples);
        float           origin[DIM];
        void            Setdx(float*, int);
        void            Setdy(float*, int);
        void            SetDirection(float*);
        double          direction[DIM];
        int             samples;
        inline float    GetRayDirection(int i){return direction[i];};
        void            PrintRayInfo();
        inline int      NumSamples(){return samples;};
        float           step;// = 254902;
};

Ray::Ray(float *o, int s)
{
    for (std::size_t i = 0; i < DIM; ++i)
        origin[i] = o[i];
    samples = s;
}

void
Ray::Setdx(float *x, int i)
{
    for (std::size_t _i = 0; _i < DIM; ++_i)
        dx[_i] = x[_i] * (2*i + 1 - WIDTH)/2.;
}

void
Ray::Setdy(float *y, int j)
{
    for (std::size_t i = 0; i < DIM; ++i)
        dy[i] = y[i] * (2*j + 1 - HEIGHT)/2.;
}

void
Ray::SetDirection(float *d)
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
    float          near, far;
    float          angle; // fov_x and fov_y
    float          W,H;
    float          position[DIM];
    float          focus[DIM];
    float          up[DIM];
    float          look[DIM];
    float          u[DIM];
    float          v[DIM];
    void           GetDelta(float*, float*, float*, float*);
    void           Pixel2Ray(Ray&, int, int);
    void           Sample(Ray&, TransferFunction&, unsigned char*, const int*, const float*, const float*, 
                          const float*, const float*);
    inline float*  GetOrigin(){return position;};
    void           ValueAt(const float*, const int*, const float*, const float*,
                           const float*, const float*, double&);
    void           InterpolateCell(const int*, const int*, const float*, const float*, 
                                   const float*, const float*, const float*,
                                   double&);
    void           Interpolate(const double, const double, const double, 
                               const double, const double, double&);
};

void
Camera::Pixel2Ray(Ray &ray, int i, int j)
{
    float dx[DIM];
    float dy[DIM];
    float direction[DIM];
    for( size_t _i = 0; _i < DIM; ++_i)
        look[_i] = focus[_i] - position[_i];
    CrossUnitNorm(look, up, u);
    //printf("ru = <%f,%f,%f>\n", u[0], u[1], u[2]);
    CrossUnitNorm(look, u, v);
    //printf("rv = <%f,%f,%f>\n", v[0], v[1], v[2]);
    // NOTE: look is normalized at this point
    for( size_t _i = 0; _i < DIM; ++_i )
    {
        dx[_i] = 2.*tan(angle/2)/W * u[_i];
        dy[_i] = 2.*tan(angle/2)/H * v[_i];
    }
    ray.Setdx(dx,i);
    ray.Setdy(dy,j);
    //printf("dx = <%f,%f,%f>\n",dx[0],dx[1],dx[2]);
    //printf("dy = <%f,%f,%f>\n",dy[0],dy[1],dy[2]);
    // Direction
    float look_norm = vtkMath::Normalize(look);
    //printf("look = <%f,%f,%f>\n", look[0], look[1], look[2]);
    //printf("look norm = %f\n", vtkMath::Normalize(look));
    direction[0] = look[0] + ((2*i+1 - W)/2)*dx[0] + ((2*j+1-H)/2)*dy[0];
    direction[1] = look[1] + ((2*i+1 - W)/2)*dx[1] + ((2*j+1-H)/2)*dy[1];
    direction[2] = look[2] + ((2*i+1 - W)/2)*dx[2] + ((2*j+1-H)/2)*dy[2];
    vtkMath::Normalize(direction);
    ray.SetDirection(direction);
    //printf("look = <%f,%f,%f>, W = %f, H = %f\n", look[0], look[1], look[2], W, H);
    //printf("Direction <%f,%f,%f>\n", direction[0], direction[1], direction[2]);
}

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
    rv.W = WIDTH;
    rv.H = HEIGHT;

    return rv;
}

void
Camera::Sample(Ray &ray, TransferFunction &tf, unsigned char *pixel, const int *dims, const float *X, const float *Y,
               const float *Z, const float *F)
{
    float sample[3] = {0,0,0};
    double offset[3] = {0,0,0};
    //float location[3] = {0,0,0};
    double val = 0;
    double dist = near;
    //unsigned char color[3] = {0,0,0};
    float opacity = 0;
    float color[3] = {0,0,0};
    for (size_t i = 0; i < DIM; ++i)
        ray.origin[i] = position[i];
    //for (size_t i = 0; i < DIM; ++i)
    //{
    //    sample[i] = ray.GetRayDirection(i) * near;
    //    location[i] = location[i] + sample[i];
    //}
    //printf("Location <%f,%f,%f>\n", location[0], location[1], location[2]);
    //ValueAt(location, dims, X, Y, Z, F, val);
    for (size_t i = 0; i < ray.samples; ++i)
    {
        for (size_t j = 0; j < DIM; ++j)
        {
            offset[j] = (double)ray.direction[j]*dist;
            sample[j] = ray.origin[j] + offset[j];
        }
            //offset[j] += ray.GetRayDirection(j)*ray.step;
        
        // Background color
        unsigned char bg[3] = {0,0,0};
        //printf("bg = <%u,%u,%u>\n", bg[0], bg[1], bg[2]);
        // Background opacity
        double bo = 0;
        //printf("ray <%f,%f,%f>\n dist %f\n", ray.direction[0], ray.direction[1],
        //        ray.direction[2], dist);
        //printf("origin <%f,%f,%f>\noffset <%f,%f,%f>\n", ray.origin[0], ray.origin[1], 
        //        ray.origin[2], offset[0], offset[1], offset[2]);
        //printf("Sample <%f,%f,%f>\n", sample[0], sample[1], sample[2]);
        //exit(0);
        ValueAt(sample, dims, X, Y, Z, F, val);
        //printf("sample[%d] <%f,%f,%f> with value %f\n", i, sample[0], sample[1], sample[2], val);
        //printf("After trans sample value %f, back color <%u,%u,%u>, backopacity %f\n",
        //        val, bg[0], bg[1], bg[2], bo);
        tf.ApplyTransferFunction(val, bg, bo);
        //printf("After trans sample value %f, back color <%u,%u,%u>, backopacity %f\n",
        //        val, bg[0], bg[1], bg[2], bo);
        bo = 1 - pow((1-bo),500./(float)ray.samples);
        //printf("Back opacity updated %f\n", bo);
        //printf("bg[0] = %u\n", bg[0]);
        color[0] = color[0] + (1-opacity)*bo*(bg[0]/255.);
        //printf("bg[0] = %u\n", bg[0]);
        //printf("color[0] = %f, opacity = %f, bo = %f, bg[0] = %f\n",
        //        color[0], opacity, bo, bg[0]);
        //printf("bo = %f\n", bo);
        //printf("%f,%f,%f,%f\n",
        //        (1-opacity), (1-opacity)*bo, (1-opacity)*bo*(bg[1]/255.),
        //        color[1] + (1-opacity)*bo*(bg[1]/255.));
        color[1] = color[1] + (1-opacity)*bo*(bg[1]/255.);
        //printf("color[1] = %f, opacity = %f, bo = %f, bg[1] = %f\n",
        //        color[1], opacity, bo, bg[1]);
        //printf("(1-opacity) %f, _*bo %f, bg[1]/255. %f\n",
        //        (1-opacity), (1-opacity)*bo, bg[1]/255.);
        color[2] = color[2] + (1-opacity)*bo*(bg[2]/255.);
        printf("Color updated <%u,%u,%u>\n", color[0], color[1], color[2]);
        opacity += (1-opacity)*bo;
        //printf("Opacity updated %f\n", opacity);
        printf("Got color <%u,%u,%u>\n", color[0], color[1], color[2]);
        dist += ray.step;
        //printf("Distance updated %f\n", dist);
        //exit(0);
        if ( i == 0) exit(0);
    }
    pixel[0] = (unsigned char)(color[0]*255);
    pixel[1] = (unsigned char)(color[1]*255);
    pixel[2] = (unsigned char)(color[2]*255);
    printf("Pixels <%u,%u,%u> from colors <%u,%u,%u>\n",
            pixel[0], pixel[1], pixel[2],
            color[0], color[1], color[2]);
}

void
Camera::ValueAt(const float *location, const int *dims, const float *X, 
                const float *Y, const float *Z, const float *F, double &val)
{
    int idx[3] = {0,0,0};
    idx[0] = BinarySearch(location[0], 0, dims[0]-1, X);
    idx[1] = BinarySearch(location[1], 0, dims[1]-1, Y);
    idx[2] = BinarySearch(location[2], 0, dims[2]-1, Z);
    InterpolateCell(idx, dims, location, X, Y, Z, F, val);
    //printf("Got value %f\n", val);
    //printf("Value at index %d = %f\n", cellIndex, F[cellIndex]);
}

void
Camera::InterpolateCell(const int *idx, const int *dims, const float *pt, const float *X, 
                        const float *Y, const float *Z, const float *F, double &val)
{
    int tmpidx[DIM];
    for(size_t i = 0; i < DIM; ++i)
        tmpidx[i] = idx[i];
    int index[8];
    double values[6];
    //printf("Orig idx %d,%d,%d\n", idx[0], idx[1], idx[2]);
    // tmpidx[1] = idx[1]
    index[0] = GetPointIndex(tmpidx, dims);  // 0,0,0
    //printf("Vertex %d,%d,%d and index %d\n", tmpidx[0], tmpidx[1], tmpidx[2], index[0]);
    tmpidx[0] = idx[0] + 1;
    index[1] = GetPointIndex(tmpidx, dims);  // 1,0,0
    //printf("Vertex %d,%d,%d and index %d\n", tmpidx[0], tmpidx[1], tmpidx[2], index[1]);
    tmpidx[0] = idx[0];
    tmpidx[2] = idx[2] + 1;
    index[2] = GetPointIndex(tmpidx, dims);  // 0,0,1
    //printf("Vertex %d,%d,%d and index %d\n", tmpidx[0], tmpidx[1], tmpidx[2], index[2]);
    tmpidx[0] = idx[0] + 1;
    tmpidx[2] = idx[2] + 1;
    index[3] = GetPointIndex(tmpidx, dims);  // 1,0,1
    //printf("Vertex %d,%d,%d and index %d\n", tmpidx[0], tmpidx[1], tmpidx[2], index[3]);
    // y + 1
    tmpidx[1] = idx[1] + 1;
    tmpidx[0] = idx[0];
    tmpidx[2] = idx[2];
    index[4] = GetPointIndex(tmpidx, dims);  // 0,1,0
    //printf("Vertex %d,%d,%d and index %d\n", tmpidx[0], tmpidx[1], tmpidx[2], index[4]);
    tmpidx[0] = idx[0] + 1;
    index[5] = GetPointIndex(tmpidx, dims);  // 1,1,0
    //printf("Vertex %d,%d,%d and index %d\n", tmpidx[0], tmpidx[1], tmpidx[2], index[5]);
    tmpidx[0] = idx[0];
    tmpidx[2] = idx[2] + 1;
    index[6] = GetPointIndex(tmpidx, dims);  // 0,1,1
    //printf("Vertex %d,%d,%d and index %d\n", tmpidx[0], tmpidx[1], tmpidx[2], index[6]);
    tmpidx[0] = idx[0] + 1;
    index[7] = GetPointIndex(tmpidx, dims);  // 1,1,1
    //printf("Vertex %d,%d,%d and index %d\n", tmpidx[0], tmpidx[1], tmpidx[2], index[7]);
    //printf("Indices %d,%d,%d,%d,%d,%d,%d\n", index[0], index[1], index[2], 
    //       index[3], index[4], index[5], index[6], index[7]);
    // Front Face
    Interpolate(pt[0], X[idx[0]], X[idx[0]+1], F[index[0]], F[index[1]],values[0]);
    Interpolate(pt[0], X[idx[0]], X[idx[0]+1], F[index[2]], F[index[3]], values[1]);
    // Connect Front Face (val 2)
    //Interpolate(pt[1], Y[idx[1]], Y[idx[1]+1], values[0], values[1], values[2]);
    Interpolate(pt[2], Z[idx[2]], Z[idx[2]+1], values[0], values[1], values[2]);
    // Back Face
    Interpolate(pt[0], X[idx[0]], X[idx[0]+1], F[index[4]], F[index[5]],values[3]);
    Interpolate(pt[0], X[idx[0]], X[idx[0]+1], F[index[6]], F[index[7]], values[4]);
    // Connect Back Face (val 5)
    //Interpolate(pt[1], Y[idx[1]], Y[idx[1]+1], values[2], values[3], values[5]);
    Interpolate(pt[2], Z[idx[2]], Z[idx[2]+1], values[3], values[4], values[5]);
    // Connect Faces
    Interpolate(pt[2], Z[idx[2]], Z[idx[2]+1], values[2], values[5], val);
    //printf("Interpolation values %f,%f,%f,%f,%f,%f,%f\n", values[0], values[1],
    //       values[2], values[3], values[4], values[5], val);
    //printf("From interpolate cell got value %f\n", val);
}

void
Camera::Interpolate(const double x, const double a, const double b,
                    const double F_A, const double F_B, double &F_X)
{
    //printf("x,a,b = %f,%f,%f\n", x, a, b);
    F_X = F_A + ( (x-a)/(b-a) )*(F_B - F_A);
}

int main()
{
    int samples = 256;
    Camera camera = SetupCamera();
    TransferFunction tf = SetupTransferFunction();
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
    //printf("Dims <%d,%d,%d>\n", dims[0], dims[1], dims[2]);
    int numPixels = dims[0]*dims[1]*dims[2];
    //unsigned char pixels[numPixels][3];

    float *X = (float*) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float*) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float*) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float*) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    // Image to create
    vtkImageData *outImage = vtkImageData::New();
    outImage->SetDimensions(camera.W, camera.H, 1);
    outImage->AllocateScalars(VTK_UNSIGNED_CHAR,3);
    unsigned char *pixels = (unsigned char*)outImage->GetScalarPointer(0,0,0);

    for (size_t i = 0; i < WIDTH; ++i )
    {
        for (size_t j = 0; j < HEIGHT; ++j)
        {
            if (i == 50 && j == 50)
            {
            Ray ray(camera.GetOrigin(), samples);// = new Ray();
            ray.step = (camera.far - camera.near) / float(ray.NumSamples()-1);
            camera.Pixel2Ray(ray, i, j);
            //printf("ray direction <%f,%f,%f>\n", ray.GetRayDirection(0), ray.GetRayDirection(1), ray.GetRayDirection(2));
            //ray.PrintRayInfo();
            unsigned char color[3];
            int index = i * j * 3;
            unsigned char pixel[3] = {pixels[index + 0], pixels[index + 1], pixels[index + 2]};
            camera.Sample(ray, tf, pixel, dims, X, Y, Z, F);
            pixels[index + 0] = pixel[0];
            pixels[index + 1] = pixel[1];
            pixels[index + 2] = pixel[2];
            printf("Got pixels <%u,%u,%u>\n", pixels[index+0], pixels[index+1], pixels[index+2]);
            printf("From <%u,%u,%u>\n", pixel[0], pixel[1], pixel[2]);
            }
        }
    }

    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(outImage);
    writer->SetFileName("output.png");
    writer->Write();
    writer->Delete();

    return 0;
}
