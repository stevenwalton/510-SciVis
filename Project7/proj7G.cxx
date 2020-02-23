
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
#include <vtkCamera.h>
#include <vtkPNGWriter.h>
#include <vtkWindowToImageFilter.h>
#include <string>

int main(int argc, char *argv[])
{
    vtkDataSetMapper *mapper = vtkDataSetMapper::New();

    vtkDataSetMapper *mapper2 = vtkDataSetMapper::New();
    
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
    vtkRenderWindow *renwin = vtkRenderWindow::New();
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    vtkDataSetReader *reader = vtkDataSetReader::New();
    vtkCutter *cut = vtkCutter::New();
    vtkPlane *plane = vtkPlane::New();
    vtkContourFilter *cf = vtkContourFilter::New();
    vtkSmartPointer<vtkPNGWriter> pngw = 
        vtkSmartPointer<vtkPNGWriter>::New();
    vtkSmartPointer<vtkImageData> img;
    vtkSmartPointer<vtkWindowToImageFilter> win2img = 
        vtkSmartPointer<vtkWindowToImageFilter>::New();
    vtkSmartPointer<vtkCamera> cam = 
        vtkSmartPointer<vtkCamera>::New();
    reader->SetFileName("noise.vtk");
    reader->Update();

    plane->SetOrigin(0,0,0);
    plane->SetNormal(0,0,1);
    cut->SetInputData(reader->GetOutput());
    cut->SetCutFunction(plane);
    cut->Update();

    cf->SetNumberOfContours(1);
    mapper->SetInputData(cut->GetOutput());
    renwin->AddRenderer(ren);
    renwin->SetSize(768,768);
    cf->SetInputData(reader->GetOutput());
    mapper2->SetInputData(cf->GetOutput());
    ren->AddActor(actor);
    ren->SetViewport(0,0,0.5,1);
    vtkRenderer *ren2 = vtkRenderer::New();
    ren2->SetViewport(0.5, 0, 1.0, 1);
    ren2->AddActor(actor2);
    renwin->AddRenderer(ren2);

    iren->SetRenderWindow(renwin);
    renwin->Render();
    for(size_t i = 0; i < 500; ++i)
    {
        double val = 1.0 + (0.01*i);
        cf->SetValue(val,val);
        cf->Update();
        ren2->GetActiveCamera()->ShallowCopy(ren->GetActiveCamera());
        iren->GetRenderWindow()->Render();
    }
    //iren->Start();
    cut->Delete();
    cf->Delete();
}


