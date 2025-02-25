/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <bitset> // For binary representation

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

class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

// Steven Functions

void
GetPointIndices(const int *index, const int *dims, int *vertex)
{
    // 2 -- 3
    // |    |
    // 0 -- 1
    // Bottom Left
    int newIndex[2] = {index[0], index[1]};
    vertex[0] = GetPointIndex(newIndex, dims);
    // Bottom Right
    newIndex[0] = index[0]+1;
    vertex[1] = GetPointIndex(newIndex, dims);
    // Upper Left
    newIndex[0] = index[0];
    newIndex[1] = index[1]+1;
    vertex[2] = GetPointIndex(newIndex, dims);
    // Upper Right
    newIndex[0] = index[0]+1;
    vertex[3] = GetPointIndex(newIndex, dims);
}

void
DetermineCase(const float *F, const float isoValue, 
              const int *vertex, int &id)
{
    // Use bitset to act in binary then convert back to an unsigned integer
    std::bitset<4> _case;
    for (size_t i = 0; i < _case.size(); ++i)
        _case[i] = F[vertex[i]] > isoValue ? 1 : 0;
    id = (int)_case.to_ulong();
}

float
InterpolateEdge(const float isoValue, const float F_a, const float F_b, 
                const float a, const float b)
{   
    return ( (isoValue - F_a)/(F_b - F_a) * (b - a)) + a;
}

void
DrawSegments(SegmentList &sl, const int x_ind, const int y_ind,
             const int *dims, const float *X, const float *Y, 
             const float *F, const int id, const float isoValue)
{
    /*
     * This is different than before
     * closer to Hank's "I've never seen and index start in the 
     * origin ;)
     * 2---1
     * |   |
     * 3---0
     */
    int vert_indices[4][2] = {{x_ind+1, y_ind},
                              {x_ind+1, y_ind+1},
                              {x_ind, y_ind+1},
                              {x_ind, y_ind}};
    float vertices[4] = {F[GetPointIndex(vert_indices[0], dims)],
                         F[GetPointIndex(vert_indices[1], dims)],
                         F[GetPointIndex(vert_indices[2], dims)],
                         F[GetPointIndex(vert_indices[3], dims)]};
    float endpoints[2][2] = {{0,0}, {0,0}};
    switch(id)
    {
        case 0:
        case 15:
            // Don't draw any lines
            break;
        case 1:
        case 14:
            endpoints[0][0] = InterpolateEdge(isoValue, vertices[3], vertices[0], X[x_ind], X[x_ind+1]);
            endpoints[0][1] = Y[y_ind];

            endpoints[1][0] = X[x_ind];
            endpoints[1][1] = InterpolateEdge(isoValue, vertices[3], vertices[2], Y[y_ind], Y[y_ind+1]);
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            break;
        case 2:
        case 13:
            endpoints[0][0] = InterpolateEdge(isoValue, vertices[3], vertices[0], X[x_ind], X[x_ind+1]);
            endpoints[0][1] = Y[y_ind];

            endpoints[1][0] = X[x_ind+1];
            endpoints[1][1] = InterpolateEdge(isoValue, vertices[1], vertices[0], Y[y_ind+1], Y[y_ind]);
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            break;
        case 3:
        case 12:
            endpoints[0][0] = X[x_ind];
            endpoints[0][1] = InterpolateEdge(isoValue, vertices[2], vertices[3], Y[y_ind+1], Y[y_ind]);

            endpoints[1][0] = X[x_ind+1];
            endpoints[1][1] = InterpolateEdge(isoValue, vertices[1], vertices[0], Y[y_ind+1], Y[y_ind]);
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            break;
        case 4:
        case 11:
            endpoints[0][0] = X[x_ind];
            endpoints[0][1] = InterpolateEdge(isoValue, vertices[3], vertices[2], Y[y_ind], Y[y_ind+1]);
                              
            endpoints[1][0] = InterpolateEdge(isoValue, vertices[1], vertices[2], X[x_ind+1], X[x_ind]);
            endpoints[1][1] = Y[y_ind+1];
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            break;
        case 5:
        case 10:
            endpoints[0][0] = InterpolateEdge(isoValue, vertices[3], vertices[0], X[x_ind], X[x_ind+1]);
            endpoints[0][1] = Y[y_ind];

            endpoints[1][0] = InterpolateEdge(isoValue, vertices[1], vertices[2], X[x_ind+1], X[x_ind]);
            endpoints[1][1] = Y[y_ind+1];
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            break;
        case 6:
            endpoints[0][0] = InterpolateEdge(isoValue, vertices[3], vertices[0], X[x_ind], X[x_ind+1]);
            endpoints[0][1] = Y[y_ind];
            endpoints[1][0] = X[x_ind+1];
            endpoints[1][1] = InterpolateEdge(isoValue, vertices[0], vertices[1], Y[y_ind], Y[y_ind+1]);
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            // Part 2
            endpoints[0][0] = X[x_ind];
            endpoints[0][1] = InterpolateEdge(isoValue, vertices[3], vertices[2], Y[y_ind], Y[y_ind+1]);
            endpoints[1][0] = InterpolateEdge(isoValue, vertices[2], vertices[1], X[x_ind], X[x_ind+1]);
            endpoints[1][1] = Y[y_ind+1];
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            break;
        case 7:
        case 8:
            endpoints[0][0] = InterpolateEdge(isoValue, vertices[1], vertices[2], X[x_ind+1], X[x_ind]);
            endpoints[0][1] = Y[y_ind+1];
            endpoints[1][0] = X[x_ind+1];
            endpoints[1][1] = InterpolateEdge(isoValue, vertices[0], vertices[1], Y[y_ind], Y[y_ind+1]);
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            break;
        case 9:
            // Other hyperbola
            endpoints[0][0] = X[x_ind];
            endpoints[0][1] = InterpolateEdge(isoValue, vertices[3], vertices[2], Y[y_ind], Y[y_ind+1]);
            endpoints[1][0] = InterpolateEdge(isoValue, vertices[3], vertices[0], X[x_ind], X[x_ind+1]);
            endpoints[1][1] = Y[y_ind];
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            // Part 2
            endpoints[0][0] = InterpolateEdge(isoValue, vertices[2], vertices[1], X[x_ind], X[x_ind+1]);
            endpoints[0][1] = Y[y_ind+1];
            endpoints[1][0] = X[x_ind+1];
            endpoints[1][1] = InterpolateEdge(isoValue, vertices[0], vertices[1], Y[y_ind], Y[y_ind+1]);
            sl.AddSegment(endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]);
            break;
        default:
            fprintf(stderr, "An error occurred. Got ID: %d", id);
            exit(1);
            break;
    }
}

void
GenIsoLines(SegmentList &sl, const int *dims, 
            const float *X, const float *Y, const float *F,
            float isoValue)
{
    int log_index[2] = {0,0};
    int vertex[4] = {0,0,0,0};
    int id;
    for (int x_ind = 0; x_ind < dims[0]-1; ++x_ind)
    {
        log_index[0] = x_ind;
        for (int y_ind = 0; y_ind < dims[1]-1; ++y_ind)
        {
            log_index[1] = y_ind;
            GetPointIndices(log_index, dims, vertex);
            DetermineCase(F, isoValue, vertex, id);
            DrawSegments(sl, x_ind, y_ind, dims, X, Y, F, id, isoValue);
        }
    }
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj5.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!
    float isoValue = 3.2; // According to the instructions.
    GenIsoLines(sl, dims, X, Y, F, isoValue);

    vtkPolyData *pd = sl.MakePolyData();
    
    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
