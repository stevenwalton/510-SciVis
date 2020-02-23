
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
//
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkContourFilter.h>

int main(int argc, char *argv[])
{
   vtkDataSetReader *reader = vtkDataSetReader::New();
   reader->SetFileName("noise.vtk");
   reader->Update();

   //
   vtkCutter *cut = vtkCutter::New();
   vtkPlane *plane = vtkPlane::New();
   plane->SetOrigin(0,0,0);
   plane->SetNormal(0,0,1);
   cut->SetInputData(reader->GetOutput());
   cut->SetCutFunction(plane);
   cut->Update();

   vtkContourFilter *cf = vtkContourFilter::New();
   cf->SetNumberOfContours(2);
   cf->SetValue(2.4,2.4);
   cf->SetValue(4,4);
   cf->SetInputData(reader->GetOutput());
   cf->Update();

   vtkDataSetMapper *mapper = vtkDataSetMapper::New();
   mapper->SetInputData(cut->GetOutput());

   vtkDataSetMapper *mapper2 = vtkDataSetMapper::New();
   mapper2->SetInputData(cf->GetOutput());
   
   vtkLookupTable *lut = vtkLookupTable::New();
   float colors[256][4];
   for (size_t i = 0; i < 255; ++i)
       lut->SetTableValue(i,i,0,255-i);
   mapper->SetLookupTable(lut);
   mapper->SetScalarRange(1,6);
   mapper2->SetLookupTable(lut);
   mapper2->SetScalarRange(1,6);
   lut->Build();

   vtkActor *actor = vtkActor::New();
   actor->SetMapper(mapper);
   vtkActor *actor2 = vtkActor::New();
   actor2->SetMapper(mapper2);

   vtkRenderer *ren = vtkRenderer::New();
   ren->AddActor(actor);
   ren->AddActor(actor2);

   vtkRenderWindow *renwin = vtkRenderWindow::New();
   renwin->AddRenderer(ren);
   renwin->SetSize(768,768);

   vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
   iren->SetRenderWindow(renwin);
   renwin->Render();
   iren->Start();

   cut->Delete();
   cf->Delete();
}


