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
#define WIDTH 1000
#define HEIGHT 1000
#define SAMPLES 1024

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
}

double
LinearInterpolation(const float a, const float b, const float x,
                    const float F_a, const float F_b)
{
    double t = (x-a)/(b-a);
    return F_a + t*(F_b - F_a);
}

int
BinarySearch(const double pt, const int lower, const int upper,
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
void CrossUnitNorm(const double *a, const double *b, double *c)
{
    vtkMath::Cross(a,b,c);
    vtkMath::Normalize(c);
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
        if (value < min || value > max) return;
        int bin = (value - min) * (256 / (max - min));
        //printf("Got bin %d\n", bin);
        RGB[0] = colors[3*bin+0];
        RGB[1] = colors[3*bin+1];
        RGB[2] = colors[3*bin+2];
        opacity = opacities[bin];
    }
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
    rv.opacities = new double[256];
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

    return rv;
}

struct Camera
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
    double          look[3];
    double          u[3];
    double          v[3];
    double          fx[3];
    double          fy[3];
    void            ViewParams();
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
    rv.angle = 30 * (PI/180.);
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;

    return rv;
}

void
Camera::ViewParams()
{
    for (size_t i = 0; i < DIM; ++i)
        look[i] = focus[i] - position[i];
    vtkMath::Normalize(look);
    CrossUnitNorm(look, up, u);
    CrossUnitNorm(look, u, v);

    for (size_t i = 0; i < DIM; ++i)
    {
        fx[i] = 2*tan(angle/2.) * u[i] / WIDTH;
        fy[i] = 2*tan(angle/2.) * v[i] / HEIGHT;
    }

}

struct Ray
{
    double      dx[3];
    double      dy[3];
    double      origin[3];
    double      direction[3];
    double      stepSize;
    void        Initialize(Camera&, const int, const int);
    void        Sample(Camera&, TransferFunction&, const int*, const float*, const float*, 
                       const float*, const float*, unsigned char*);
    void        GetValue(const double*, const int*, const float*, const float*,
                         const float*, const float*, double&);
    void        GenIndices(const double*, const int*, const float*, const float*, 
                           const float*, int*);
    void        LERPCell(const double*, const int*, const float*, const float*,
                        const float*, const float*, const int*, double &);
};

void
Ray::Initialize(Camera &camera, const int pixel_x, const int pixel_y)
{
    for (size_t i = 0; i < DIM; ++i)
    {
        dx[i] = camera.fx[i] * (2*pixel_x+1 - WIDTH)/2.;
        dy[i] = camera.fy[i] * (2*pixel_y+1 - HEIGHT)/2.;
        direction[i] = camera.look[i] + dx[i] + dy[i];
        origin[i] = camera.position[i];
    }
    vtkMath::Normalize(dx);
    vtkMath::Normalize(dy);
    vtkMath::Normalize(direction);
}

void
Ray::Sample(Camera &camera, TransferFunction &tf, const int *dims, const float *X, const float *Y,
            const float *Z, const float *F, unsigned char *pixel)
{
    for (size_t i = 0; i < DIM; ++i)
        pixel[i] = 0;
    double dist = camera.near;
    double sample[3] = {0,0,0};
    double color[3] = {0,0,0};
    double opacity = 0;
    for (size_t i = 0; i < SAMPLES; ++i)
    {
        double offset[3];
        double sampleValue = 0;
        double bg_opacity = 0;
        unsigned char bg_color[3] = {0,0,0};
        for (size_t j = 0; j < DIM; ++j)
        {
            offset[j] = direction[j]*dist;
            sample[j] = origin[j] + offset[j];
        }
        //printf("Origin <%f,%f,%f>\nOffset <%f,%f,%f>\n",
        //        origin[0], origin[1], origin[2],
        //        offset[0], offset[1], offset[2]);
        //printf("Sampling at <%f,%f,%f>\n", sample[0], sample[1], sample[2]);
        GetValue(sample, dims, X, Y, Z, F, sampleValue);
        //printf("Got sample value %f\n", sampleValue);
        tf.ApplyTransferFunction(sampleValue, bg_color, bg_opacity);
        bg_opacity = 1 - pow( (1-bg_opacity), 500. / (double)SAMPLES);

        color[0] = color[0] + (1-opacity) * bg_opacity * (bg_color[0] / 255.);
        color[1] = color[1] + (1-opacity) * bg_opacity * (bg_color[1] / 255.);
        //printf("Color[1] = %u\n", color[1]);
        color[2] = color[2] + (1-opacity) * bg_opacity * (bg_color[2] / 255.);
        opacity = opacity + (1-opacity)*bg_opacity;
        dist += stepSize;
        //printf("Sample[%d] was mapped by transfer function to %u,%u,%u\n opacity is now %f\n bg_opacity is %f\n", i, 
        //        bg_color[0], bg_color[1], bg_color[2], opacity, bg_opacity);
        //printf("Got color <%f,%f,%f>\n", color[0], color[1], color[2]);
        //printf("Opacity %f\n", opacity);
        //if (i == 16) exit(0);
        //if (i == 21) exit(0);
    }
    pixel[0] = (unsigned char)(color[0] * 255);
    pixel[1] = (unsigned char)(color[1] * 255);
    pixel[2] = (unsigned char)(color[2] * 255);
    //printf("Got colors <%u,%u,%u>\n", pixel[0], pixel[1], pixel[2]);
}

void
Ray::GetValue(const double *pt, const int *dims, 
              const float *X, const float *Y, const float *Z, const float *F,
              double &val)
{
    int idx[3];
    GenIndices(pt, dims, X, Y, Z, idx);
    //printf("Got vertex <%d,%d,%d>\n", idx[0], idx[1], idx[2]);
    LERPCell(pt, dims, X, Y, Z, F, idx, val);
    //printf("Value = %f\n",val);
}

void 
Ray::GenIndices(const double *pt, const int *dims, const float *X, const float *Y, 
                const float *Z, int *vertex)
{
    vertex[0] = BinarySearch(pt[0], 0, dims[0]-1, X);
    vertex[1] = BinarySearch(pt[1], 0, dims[1]-1, Y);
    vertex[2] = BinarySearch(pt[2], 0, dims[2]-1, Z);
}

void
Ray::LERPCell(const double *pt, const int *dims, const float *X, const float *Y,
              const float *Z, const float *F, const int *idx, double &val)
{
    int tmpidx[3];
    for (size_t i = 0; i < DIM; ++i)
        tmpidx[i] = idx[i];
    int index[8];
    double values[6];

    // Front Face y = 0
    index[0] = GetPointIndex(tmpidx, dims); // 0,0,0
    tmpidx[0] = idx[0] + 1;
    index[1] = GetPointIndex(tmpidx, dims); // 1,0,0
    tmpidx[0] = idx[0];
    tmpidx[2] = idx[2] + 1;
    index[2] = GetPointIndex(tmpidx, dims); // 0,0,1
    tmpidx[0] = idx[0] + 1;
    index[3] = GetPointIndex(tmpidx, dims); // 1,0,1
    // Back Face y = 1
    tmpidx[0] = idx[0];
    tmpidx[1] = idx[1] + 1;
    tmpidx[2] = idx[2];
    index[4] = GetPointIndex(tmpidx, dims); // 0,1,0
    tmpidx[0] = idx[0] + 1;
    index[5] = GetPointIndex(tmpidx, dims); // 1,1,0
    tmpidx[0] = idx[0];
    tmpidx[2] = idx[2] + 1;
    index[6] = GetPointIndex(tmpidx, dims); // 0,1,1
    tmpidx[0] = idx[0] + 1;
    //tmpidx[2] = idx[2] + 1;
    index[7] = GetPointIndex(tmpidx, dims); // 1,1,1


    // Interpolate
    values[0] = LinearInterpolation(X[idx[0]], X[idx[0]+1], pt[0], F[index[0]], F[index[1]]);
    values[1] = LinearInterpolation(X[idx[0]], X[idx[0]+1], pt[0], F[index[2]], F[index[3]]);
    // combine
    values[2] = LinearInterpolation(Z[idx[2]], Z[idx[2]+1], pt[2], values[0], values[1]);
    // Interpolate
    values[3] = LinearInterpolation(X[idx[0]], X[idx[0]+1], pt[0], F[index[4]], F[index[5]]);
    values[4] = LinearInterpolation(X[idx[0]], X[idx[0]+1], pt[0], F[index[6]], F[index[7]]);
    // combine
    values[5] = LinearInterpolation(Z[idx[2]], Z[idx[2]+1], pt[2], values[3], values[4]);
    // Full combine
    val = LinearInterpolation(Y[idx[1]], Y[idx[1]+1], pt[1], values[2], values[5]);
}


int main()
{
    int dims[3];

    // Read in the DataSet, which is a Rectilinear/Uniform Grid
    vtkDataSetReader *reader = vtkDataSetReader::New();
    //reader->SetFileName("astro64.vtk");
    reader->SetFileName("astro512.vtk");
    reader->Update();
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *)reader->GetOutput();
    rgrid->GetDimensions(dims);

    vtkImageData *outImage = vtkImageData::New();
    outImage->SetDimensions(WIDTH, HEIGHT, 1);
    outImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    unsigned char *outImgData = (unsigned char*)outImage->GetScalarPointer(0,0,0);

    float *X = (float *)rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *)rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *)rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *)rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    Camera camera = SetupCamera();
    camera.ViewParams();
    TransferFunction tf = SetupTransferFunction();

    Ray ray;
    ray.stepSize = (camera.far - camera.near) / (SAMPLES - 1);
    for (size_t i = 0; i < WIDTH; ++i)
    {
        for (size_t j = 0; j < HEIGHT; ++j)
        {
            //if (i == 50 && j == 50)
            //if (i == 19 && j == 76)
            //{
            ray.Initialize(camera,i,j);
            //printf("ray <%f,%f,%f>\n", ray.direction[0], ray.direction[1], ray.direction[2]);
            unsigned char pixelColor[3] = {0,0,0};
            ray.Sample(camera, tf, dims, X, Y, Z, F, pixelColor);
            //exit(0);
            int index = (j*WIDTH + i)*3;
            outImgData[index + 0] = pixelColor[0];
            outImgData[index + 1] = pixelColor[1];
            outImgData[index + 2] = pixelColor[2];
            //printf("Imagedata[%d] = %u\n", index, outImgData[index]);
            //printf("Imagedata[%d] = %u\n", index+1, outImgData[index+1]);
            //printf("Imagedata[%d] = %u\n", index+2, outImgData[index+2]);
            //printf("Got out img data <%u,%u,%u>\nFrom <%u,%u,%u>\n",
            //        outImgData[index], outImgData[index+1], outImgData[index+2],
            //        pixelColor[0], pixelColor[1], pixelColor[2]);
            //}
            //if (i == 19 && j == 76) exit(0);
        }
    }
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(outImage);
    writer->SetFileName("output.png");
    writer->Write();
    writer->Delete();

}
