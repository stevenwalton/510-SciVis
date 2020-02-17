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
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <bitset>
#include <iostream>

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
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
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
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
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
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
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
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
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
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    // idx[0] = pointId%dims[0];
    // idx[1] = pointId/dims[0];
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
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}


class TriangleList
{
   public:
                   TriangleList() { maxTriangles = 1000000; triangleIdx = 0; pts = new float[9*maxTriangles]; };
     virtual      ~TriangleList() { delete [] pts; };

     void          AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxTriangles;
     int           triangleIdx;
};

void
TriangleList::AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3)
{
    pts[9*triangleIdx+0] = X1;
    pts[9*triangleIdx+1] = Y1;
    pts[9*triangleIdx+2] = Z1;
    pts[9*triangleIdx+3] = X2;
    pts[9*triangleIdx+4] = Y2;
    pts[9*triangleIdx+5] = Z2;
    pts[9*triangleIdx+6] = X3;
    pts[9*triangleIdx+7] = Y3;
    pts[9*triangleIdx+8] = Z3;
    triangleIdx++;
}

vtkPolyData *
TriangleList::MakePolyData(void)
{
    int ntriangles = triangleIdx;
    int numPoints = 3*(ntriangles);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *tris = vtkCellArray::New();
    tris->EstimateSize(numPoints,4);
    for (int i = 0 ; i < ntriangles ; i++)
    {
        double pt[3];
        pt[0] = pts[9*i];
        pt[1] = pts[9*i+1];
        pt[2] = pts[9*i+2];
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[9*i+3];
        pt[1] = pts[9*i+4];
        pt[2] = pts[9*i+5];
        vtk_pts->SetPoint(ptIdx+1, pt);
        pt[0] = pts[9*i+6];
        pt[1] = pts[9*i+7];
        pt[2] = pts[9*i+8];
        vtk_pts->SetPoint(ptIdx+2, pt);
        vtkIdType ids[3] = { ptIdx, ptIdx+1, ptIdx+2 };
        tris->InsertNextCell(3, ids);
        ptIdx += 3;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetPolys(tris);
    tris->Delete();
    vtk_pts->Delete();

    return pd;
}

class Tetrahedron
{
  public:
    float X[4];
    float Y[4];
    float Z[4];
    float F[4];
    void PrintSelf()
    {
        for (int i = 0 ; i < 4 ; i++)
            printf("\tV%d: (%f, %f, %f) = %f\n", i, X[i], Y[i], Z[i], F[i]);
    };
};

float
InterpolateEdge(const float isoValue, const float F_a, const float F_b, 
                const float a, const float b)
{   
    return ( (isoValue - F_a)/(F_b - F_a) * (b - a)) + a;
}

void
MakeVectors(float pts[4][3], float vec[4][3])
{
    // X
    vec[0][0] = pts[1][0] - pts[0][0];
    vec[1][0] = pts[2][0] - pts[1][0];
    vec[2][0] = pts[3][0] - pts[2][0];
    vec[3][0] = pts[0][0] - pts[3][0];
    // Y
    vec[0][1] = pts[1][1] - pts[0][1];
    vec[1][1] = pts[2][1] - pts[1][1];
    vec[2][1] = pts[3][1] - pts[2][1];
    vec[3][1] = pts[0][1] - pts[3][1];
    // Z
    vec[0][2] = pts[1][2] - pts[0][2];
    vec[1][2] = pts[2][2] - pts[1][2];
    vec[2][2] = pts[3][2] - pts[2][2];
    vec[3][2] = pts[0][2] - pts[3][2];
}

void
Hypotnuse(float a[3], float b[3], float c[3])
{
    c[0] = a[0]*a[0] + b[0]*b[0];
    c[1] = a[1]*a[1] + b[1]*b[1];
    c[2] = a[2]*a[2] + b[2]*b[2];
}

void 
IsosurfaceTet(Tetrahedron &tet, TriangleList &tl, float isoval)
{
    float eps = 0.01;
    int sum = 0;
    std::bitset<4> _case;
    // Determine cases
    for (size_t i = 0; i < 4; ++i)
    {
        sum += tet.F[i] > isoval ? 1 : 0; 
        _case[i] = tet.F[i] > isoval ? 1 : 0;
    }
    float pts[4][3]; // We have potentially 4 points to use
    size_t point = 0;
    // Return if we aren't drawing the tet
    switch (sum)
    {
        // All F > isoval || All F < isoval
        case 0:
        case 4:
            return;
        // 1 F > isoval || 1 F < isoval
        case 1:
        case 3:
        // 2 F > isoval && 2 F < isoval
        // Special case where we have to make 2 triangles
        case 2:
        {
            point = 0;
            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = i+1; j < 4; ++j)
                {
                    if (( _case[i] == 1 && _case[j] == 0) ||
                        ( _case[i] == 0 && _case[j] == 1))
                    {
                        pts[point][0] = InterpolateEdge(isoval, tet.F[i], tet.F[j],
                                                        tet.X[i], tet.X[j]);
                        pts[point][1] = InterpolateEdge(isoval, tet.F[i], tet.F[j],
                                                        tet.Y[i], tet.Y[j]);
                        pts[point][2] = InterpolateEdge(isoval, tet.F[i], tet.F[j],
                                                        tet.Z[i], tet.Z[j]);
                        if (point == 2)
                        {
                            tl.AddTriangle(pts[0][0], pts[0][1], pts[0][2],
                                           pts[1][0], pts[1][1], pts[1][2],
                                           pts[2][0], pts[2][1], pts[2][2]);
                            // Return if case 1 or 3
                            if (sum != 2) return;
                        }
                        else if ( point == 3)
                        {
                            float myVec[4][3];
                            float abc[3][3];
                            MakeVectors(pts, myVec);

                            float hyp[3];
                            Hypotnuse(myVec[0], myVec[1], hyp);
                            if (abs(hyp[0] - myVec[2][0]) < eps &&
                                abs(hyp[1] - myVec[2][1]) < eps &&
                                abs(hyp[2] - myVec[2][2]) < eps)
                            {
                                printf("option 1\n");
                                tl.AddTriangle(pts[0][0], pts[0][1], pts[0][2],
                                               pts[2][0], pts[2][1], pts[2][2],
                                               pts[3][0], pts[3][1], pts[3][2]);
                            }
                            else if (abs(hyp[0] - myVec[0][0]) < eps &&
                                     abs(hyp[1] - myVec[0][1]) < eps &&
                                     abs(hyp[2] - myVec[0][2]) < eps)
                            {
                                printf("Option 2\n");
                                tl.AddTriangle(pts[0][0], pts[0][1], pts[0][2],
                                               pts[1][0], pts[1][1], pts[1][2],
                                               pts[3][0], pts[3][1], pts[3][2]);
                            }
                            else if (abs(hyp[0] - myVec[1][0]) < eps &&
                                     abs(hyp[1] - myVec[1][1]) < eps &&
                                     abs(hyp[2] - myVec[1][2]) < eps)
                            {
                                printf("Option 3\n");
                                tl.AddTriangle(pts[1][0], pts[1][1], pts[1][2],
                                               pts[3][0], pts[3][1], pts[3][2],
                                               pts[2][0], pts[2][1], pts[2][2]);
                            }
                            else
                            {
                                printf("We have a fourth option?\n");
                                //printf("Vecsum (%f,%f,%f)\n", vecsum[0], vecsum[1], vecsum[2]);
                                printf("Have points, in order \
                                        \n(%2.2f,%f,%f), (%f,%f,%f), (%f,%f,%f), (%f,%f,%f)\n",
                                        pts[0][0], pts[0][1], pts[0][2],
                                        pts[1][0], pts[1][1], pts[1][2],
                                        pts[2][0], pts[2][1], pts[2][2],
                                        pts[3][0], pts[3][1], pts[3][2]);
                                printf("a: (%f,%f,%f)\n", myVec[0][0], myVec[0][1], myVec[0][2]);
                                printf("b: (%f,%f,%f)\n", myVec[1][0], myVec[1][1], myVec[1][2]);
                                printf("c: (%f,%f,%f)\n", myVec[2][0], myVec[2][1], myVec[2][2]);
                                printf("Hypotnuse: %f %f %f\n", hyp[0], hyp[1], hyp[2]);
                                abort();
                            }

                            return;
                        }
                        else if ( point > 4)
                        {
                            fprintf(stderr, "We got %lu points!\n", point);
                            abort();
                        }
                        point++;
                    }

                }
            }
            return;
            break;
        }
        default:
            fprintf(stderr, "An error occurred and we got a tetrahedron that \
                    cannot exist\n");
            abort();
            break;
    }
}

/*
void
DetermineCase(Tetrahedron &tet, float isoval, float F, int &id)
{
    std::bitset<4> _case;
    for (size_t i = 0; i < _case.size(); ++i)
        _case[i] = F[vertex[i]] > isoval ? 1 : 0;
    id = (int)_case.to_ulong();
}
*/

int main()
{
    int  i, j;

/* If you want to try on one tetrahedron.
    Tetrahedron t;
    t.X[0] = 0;
    t.Y[0] = 0;
    t.Z[0] = 0;
    t.F[0] = 0.6;
    t.X[1] = 1;
    t.Y[1] = 0;
    t.Z[1] = 0;
    t.F[1] = 0;
    t.X[2] = 0;
    t.Y[2] = 1;
    t.Z[2] = 0;
    t.F[2] = 0;
    t.X[3] = 0.5;
    t.Y[3] = 0.5;
    t.Z[3] = 1.0;
    t.F[3] = 1;
    TriangleList tl;
    IsosurfaceTet(t, tl, 0.5);
*/

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6.vtk");
    rdr->Update();
    if (rdr->GetOutput() == NULL || rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Could not find input file." << endl;
        exit(EXIT_FAILURE);
    }

    vtkUnstructuredGrid *ugrid = (vtkUnstructuredGrid *) rdr->GetOutput();
    float *pts = (float *) ugrid->GetPoints()->GetVoidPointer(0);
    float *F = (float *) ugrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    vtkCellArray *cell_array = ugrid->GetCells();

    TriangleList tl;
    cell_array->InitTraversal();
    int ncells = cell_array->GetNumberOfCells();
    cerr << "Number of cells to tetrahedralize is " << ncells << endl;
    int cnt = 0;
    float isoval = 12.2;
    for (int i = 0 ; i < ncells ; i++)
    {
        vtkIdType npts;
        vtkIdType *ids;
        cell_array->GetNextCell(npts, ids);
        if (npts == 4)
        {
            Tetrahedron tet;
            for (int j = 0 ; j < 4 ; j++)
            {
                // This data set is in a huge bounding box.  Normalize as we go.
                tet.X[j] = (pts[3*ids[j]]+3e+7)/6e+7;
                tet.Y[j] = (pts[3*ids[j]+1]+3e+7)/6e+7;
                tet.Z[j] = (pts[3*ids[j]+2]+3e+7)/6e+7;
                tet.F[j] = F[ids[j]];
            }
            IsosurfaceTet(tet, tl, isoval);
        }
        else
        {
            cerr << "Input was non-tetrahedron!!  Ignoring..." << endl;
            cerr << "Type is " << npts << endl;
            cerr << "Cell is " << i << endl;
            cnt++;
            continue;
        }
    }

    vtkPolyData *pd = tl.MakePolyData();

/*
    //This can be useful for debugging
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("proj6_out.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkCleanPolyData *cpd = vtkCleanPolyData::New();
    cpd->SetInputData(pd);
    cpd->SetAbsoluteTolerance(0);
    cpd->PointMergingOn();
    cpd->Update();
    vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
    pdn->SetInputData(cpd->GetOutput());
    //pdn->SetInputData(pd);
    pdn->Update();

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pdn->GetOutput());
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

    ren1->GetActiveCamera()->SetFocalPoint(0.5, 0.5, 0.5);
    ren1->GetActiveCamera()->SetPosition(0,0,-2);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(0.01, 4);
    //ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
