#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

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
    if (pt > arr[index])
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

float
LinearInterpolation(const float a, const float b, const float x,
                    const float F_a, const float F_b)
{
    float t = (x-a)/(b-a);
    return F_a + t*(F_b - F_a);
}


// ****************************************************************************
//  Function: EvaluateVectorFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//        The third value in the field is the z-component, which is 0.
//
//     rv (output): the interpolated field value. (0,0) if the location is out of bounds.
//
// ****************************************************************************

void EvaluateVectorFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F, float *rv)
{
    // IMPLEMENT ME!

    if (pt[0] < X[0] || pt[0] > X[dims[0]-1] || 
        pt[1] < Y[0] || pt[1] > Y[dims[1]-1])
    {
        rv[0] = 0; // setting the x-component of the velocity
        rv[1] = 0; // setting the y-component of the velocity
    }
    else
    {
        int idx[2];
        idx[0] = BinarySearch(pt[0], 0, dims[0]-1, X);
        idx[1] = BinarySearch(pt[1], 0, dims[1]-1, Y);
        int index = 2*idx[1]*dims[0] + 2*idx[0];
        int upper_index = 2*(idx[1]+1)*dims[0] + 2*idx[0];
        // Bottom x-line interpolation
        float l2rx = LinearInterpolation(X[idx[0]], X[idx[0]+1], pt[0],
                                         F[index], F[index+2]);
        float l2ry = LinearInterpolation(X[idx[0]], X[idx[0]+1], pt[0],
                                         F[index+1], F[index+3]);
        // Top x-line interpolation
        float t_l2rx = LinearInterpolation(X[idx[0]], X[idx[0]+1], pt[0],
                                           F[upper_index], F[upper_index+2]);
        float t_l2ry = LinearInterpolation(X[idx[0]], X[idx[0]+1], pt[0],
                                           F[upper_index+1], F[upper_index+3]);
        // Y-line interpolation
        rv[0] = LinearInterpolation(Y[idx[1]], Y[idx[1]+1], pt[1],
                                    l2rx, t_l2rx);
        rv[1] = LinearInterpolation(Y[idx[1]], Y[idx[1]+1], pt[1],
                                          l2ry, t_l2ry);
    }
        
}

// ****************************************************************************
//  Function: AdvectWithEulerStep
//
//  Arguments:
//     pt: the seed location (two-dimensions)
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//     h: The size of the Euler step
//     nsteps: The number of Euler steps to take
//     output_locations (output): An array of size 2*(nsteps+1).  It's first entry
//        should be the seed location.  The second entry should be the result
//        of the first advection step.  The final entry should be the result
//        of the final advection step.
//     speeds (output): An array of size (nsteps+1).  It's first entry should be the
//        speed at the seed location.  It's final entry should be the speed
//        at the location of the result of the final advection step.
//        Recall that speed is the magnitude of the velocity.
//
// ****************************************************************************

void
AdvectWithEulerStep(const float *pt, const int *dims, const float *X, 
                    const float *Y, const float *F, 
                    float h, int nsteps, float *output_locations, float *speeds)
{
    // Initialization
    output_locations[0] = pt[0];
    output_locations[1] = pt[1];
    float npt[2] = {pt[0], pt[1]};
    float speed[2] = {0,0};
    EvaluateVectorFieldAtLocation(pt, dims, X, Y, F, speed);
    speeds[0] = sqrt((speed[0]*speed[0]) + (speed[1]*speed[1]));
    for (int i = 1; i <= nsteps; ++i)
    {
        // Get speed at starting position
        EvaluateVectorFieldAtLocation(npt, dims, X, Y, F, speed);
        // Update speeds array with magnitude of speed
        speeds[i] = sqrt((speed[0]*speed[0]) + (speed[1]*speed[1]));
        // Get new location
        npt[0] = npt[0] + h*speed[0];
        npt[1] = npt[1] + h*speed[1];
        // Update output_locations array
        output_locations[i*2]   = npt[0];
        output_locations[i*2+1] = npt[1];
    }
}

void
AdvectWithRK4Step(const float *pt, const int *dims, const float *X, 
                  const float *Y, const float *F, 
                  float h, int nsteps, float *output_locations, float *speeds)
{

    output_locations[0] = pt[0];
    output_locations[1] = pt[1];
    float npt[2] = {pt[0], pt[1]};  // Loop positions
    float tpt[2] = {pt[0], pt[1]};  // rk positions

    // Velocity (k's)
    float velocity[4][2] = {{0,0},{0,0},{0,0},{0,0}};

    // First step: i=0
    EvaluateVectorFieldAtLocation(pt, dims, X, Y, F, velocity[0]);
    float v[2] = {velocity[0][0], velocity[0][1]};
    speeds[0] = sqrt(v[0]*v[0] + v[1]*v[1]);
    for (int i = 1; i <= nsteps; ++i) // f_n( f_{n-1} )
    {
        // K0
        tpt[0] = npt[0];
        tpt[1] = npt[1];
        EvaluateVectorFieldAtLocation(tpt, dims, X, Y, F, velocity[0]);

        // K1
        tpt[0] = npt[0] + h/2.f;
        tpt[1] = npt[1] + velocity[0][1]/2.f;
        EvaluateVectorFieldAtLocation(tpt, dims, X, Y, F, velocity[1]);

        // K2
        tpt[1] = npt[1] + velocity[1][1]/2.f;
        EvaluateVectorFieldAtLocation(tpt, dims, X, Y, F, velocity[2]);

        // K3
        tpt[0] = npt[0] + h;
        tpt[1] = npt[1] + velocity[2][1];
        EvaluateVectorFieldAtLocation(tpt, dims, X, Y, F, velocity[3]);

        // Update velocity
        v[0] = (velocity[0][0] 
                + 2.f*velocity[1][0] 
                + 2.f*velocity[2][0] 
                + velocity[3][0])
               / 6.f;
        v[1] = (velocity[0][1] 
                + 2.f*velocity[1][1] 
                + 2.f*velocity[2][1] 
                + velocity[3][1])
               / 6.f;

        // Update speed (magnitude of velocity)
        speeds[i] = sqrt(v[0]*v[0] + v[1]*v[1]);
        // Update starting points
        npt[0] = npt[0] + h*v[0];
        npt[1] = npt[1] + h*v[1];
        // Update output locations
        output_locations[i*2] = npt[0];
        output_locations[i*2+1] = npt[1];
        
        // Working Euler
        /*
        v[0] = velocity[0][0];
        v[1] = velocity[0][1];
        speeds[i] = sqrt(v[0]*v[0] + v[1]*v[1]);
        npt[0] = npt[0] + h*v[0];
        npt[1] = npt[1] + h*v[1];
        output_locations[i*2] = npt[0];
        output_locations[i*2+1] = npt[1];
        */
    }
}


// ****************************************************************************
//  Function: CalculateArcLength
//
//  Arguments:
//     locations: an array of 2D locations.
//     nlocations: the number of locations in the array "locations".
//
//  Returns: the arc length, meaning the distance between each successive
//           pair of points
//
// ****************************************************************************

float
CalculateArcLength(const float *output_locations, int nlocations)
{
    float length = 0;
    float dx = 0;
    float dy = 0;
    bool print = false;
    if (output_locations[0] > 8.623 && output_locations[0] < 8.624)
        print = true;
    for (int i = 1; i < nlocations; ++i)
    {
        dx = output_locations[i*2] - output_locations[i*2-2];
        dy = output_locations[i*2+1] - output_locations[i*2-1];
        // Get Magnitude
        length += sqrt(dx*dx + dy*dy);
    }
    return length;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// VTK files are only 3D, so the vector data is all of the form (X,Y,0).
// Remove the 0's since it is counter-intuitive for students who are 
// thinking of this as 2D data.
float *
Convert3DVectorDataTo2DVectorData(const int *dims, const float *F)
{
    float *rv = new float[dims[0]*dims[1]*2];
    int index3D = 0;
    int index2D = 0;
    for (int i = 0 ; i < dims[0] ; i++)
       for (int j = 0 ; j < dims[1] ; j++)
       {
           rv[index2D]   = F[index3D];
           rv[index2D+1] = F[index3D+1];
           index2D += 2;
           index3D += 3;
       }

    return rv;
}

vtkPolyData *
CreateVTKPolyData(int nseeds, int nsteps, float **output_locations, float **speeds)
{
    int numPoints = nseeds*(nsteps+1);
    vtkPoints *pts = vtkPoints::New();
    pts->SetNumberOfPoints(numPoints);
    vtkFloatArray *var = vtkFloatArray::New();
    var->SetName("speed");
    var->SetNumberOfTuples(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nseeds ; i++)
    {
        for (int j = 0 ; j < nsteps+1 ; j++)
        {
            double pt[3];
            pt[0] = output_locations[i][2*j];
            pt[1] = output_locations[i][2*j+1];
            pt[2] = 0.;
            pts->SetPoint(ptIdx, pt);
            var->SetTuple1(ptIdx, speeds[i][j]);
            if (j > 0)
            {
                vtkIdType ids[2] = { ptIdx-1, ptIdx };
                lines->InsertNextCell(2, ids);
            }
            ptIdx++;
        }
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(pts);
    pd->GetPointData()->AddArray(var);
    pd->GetPointData()->SetActiveScalars("speed");
    pd->SetLines(lines);
    lines->Delete();
    var->Delete();
    pts->Delete();

    return pd;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj4_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F_3D = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);
    float *F = Convert3DVectorDataTo2DVectorData(dims, F_3D);
    
    
    const int npts = 10;
    float pt[npts][3] =
         {
            {10.1119, 1.22062, 0},
            {8.62376, 13.3839, 0},
            {1.55026, 1.26123, 0},
            {6.9736, 0.653565, 0},
            {2, 2.74117, 0},
            {8.93699, 10.4111, 0},
            {6.08791, -0.533753, 0},
            {10.0543, 1.38024, 0},
            {3.84128, -0.768977, 0},
            {6.66757, 6.0259, 0},
         };


    for (i = 0 ; i < npts ; i++)
    {
       float vec[2];
       EvaluateVectorFieldAtLocation(pt[i], dims, X, Y, F, vec);
       cerr << "Velocity at (" << pt[i][0] <<", "<<pt[i][1] << ") is (" << vec[0] << ", " << vec[1] << ")" << endl;
    }

    //float h = 0.01;
    float h = 5.;
    const int nsteps = 5000;
    float **euler_output_locations = new float*[2*(npts+1)];
    float **rk4_output_locations   = new float*[2*(npts+1)];
    float **euler_speeds = new float*[npts+1];
    float **rk4_speeds   = new float*[npts+1];
    for (i = 0 ; i < npts ; i++)
    {
       euler_output_locations[i] = new float[(nsteps+1)*2];
       rk4_output_locations[i]   = new float[(nsteps+1)*2];
       euler_speeds[i] = new float[nsteps];
       rk4_speeds[i]   = new float[nsteps];
       AdvectWithEulerStep(pt[i], dims, X, Y, F, h, nsteps, euler_output_locations[i], euler_speeds[i]);
       AdvectWithRK4Step(pt[i], dims, X, Y, F, h, nsteps, rk4_output_locations[i], rk4_speeds[i]);
       float euler_length = CalculateArcLength(euler_output_locations[i], nsteps+1);
       float rk4_length = CalculateArcLength(rk4_output_locations[i], nsteps+1);
       cerr << "Arc length for (" << pt[i][0] << ", " << pt[i][1] << ") is " 
            << "Euler: " << euler_length <<  " RK4: " << rk4_length << endl;
    }

    //vtkPolyData *pd = CreateVTKPolyData(npts, nsteps, output_locations, speeds);
    vtkPolyData *euler_pd = CreateVTKPolyData(npts, nsteps, euler_output_locations, euler_speeds);
    vtkPolyData *rk4_pd = CreateVTKPolyData(npts, nsteps, rk4_output_locations, rk4_speeds);

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInput(pd);
    writer->Write();
 */

    // Render 1 - Euler
    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(euler_pd);
    win1Mapper->SetScalarRange(0, 0.15);
  
    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);
  
    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();
  
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
  
    ren1->GetActiveCamera()->SetFocalPoint(5,5,0);
    ren1->GetActiveCamera()->SetPosition(5,5,30);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);
    // Placement in Render
    ren1->SetViewport(0,0,0.5,1);

    // Render 2 - RK4
    vtkSmartPointer<vtkDataSetMapper> win2Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win2Mapper->SetInputData(rk4_pd);
    win2Mapper->SetScalarRange(0, 0.15);
  
    vtkSmartPointer<vtkActor> win2Actor =
      vtkSmartPointer<vtkActor>::New();
    win2Actor->SetMapper(win2Mapper);
  
    vtkSmartPointer<vtkRenderer> ren2 =
      vtkSmartPointer<vtkRenderer>::New();

    ren2->GetActiveCamera()->SetFocalPoint(5,5,0);
    ren2->GetActiveCamera()->SetPosition(5,5,30);
    ren2->GetActiveCamera()->SetViewUp(0,1,0);
    ren2->GetActiveCamera()->SetClippingRange(20, 120);
    ren2->GetActiveCamera()->SetDistance(30);
    // 
    ren2->AddActor(win2Actor);
    ren2->SetBackground(0.0, 0.0, 0.0);
    // Placement in Render
    ren2->SetViewport(0.5, 0, 1.0, 1);

  
    // Renderer
    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();

    renWin->AddRenderer(ren1);
    renWin->AddRenderer(ren2);
    renWin->SetWindowName("Euler vs RK4");
  
    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    renWin->SetSize(800, 800);
  

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    delete [] F;
    //pd->Delete();
    euler_pd->Delete();
    rk4_pd->Delete();
}
