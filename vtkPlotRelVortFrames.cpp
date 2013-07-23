// Pete Bosler
// March 28, 2012
// Animate .vtk files 
// filename : vtkPlotRelVortFrames.cpp
// Last updated : 3-28-12
#include "vtkActor.h"
#include "vtkAxesActor.h"
#include "vtkCamera.h"
#include "vtkFollower.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkProperty2D.h"
#include "vtkCaptionActor2D.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkScalarBarActor.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkLookupTable.h"
#include "vtkWindowToImageFilter.h"
#include "vtkJPEGWriter.h"
#include <iostream>
#include <cmath>
#include <string.h>
#include <sstream>

using namespace std;

int main( int argc, const char *argv[]){
	if ( argc != 3){
		cout << "Usage : "<< argv[0] << " <jobtitle | frame total>\n";
		return 1;
	};
	char fileName[128];
	char counterString[4];
	int j;
	int counter = 0;
	int finalFrame;
	char jpgFileName[256];
	char jpgFileRoot[128];
// USER DATA (alter this section): 
//TO DO: User input as arguments or read from data files
	float dt = 0.01;  // time step size
	double min = -9.0;	// minimum tracer value
	double max = 9.0;   // maximum tracer value
	char scalarsDataName[16]="relVortPanel";	// vtk data name (must match .vtk file)
	bool edgesOn = true;	// Plot panel edges (y/n)
	char scalarName[64] = "Rel. Vort.";	// Title displayed in plot


// Set up input filenames
	j=sprintf(counterString,"%04d",counter);
	strcpy(fileName,argv[1]);
	strcat(fileName,counterString);
	strcat(fileName,".vtk");
// Set up output filenames	
	strcpy(jpgFileName,"jpgOut/relVortMovie_");
	strcpy(jpgFileRoot,jpgFileName);
	strcat(jpgFileName,counterString);
	strcat(jpgFileName,".jpg");
	cout << " Plotting files "<< fileName << " to ";
	istringstream ss(argv[2]);
	if (!(ss>>finalFrame)) cerr << "Invalid number "<<argv[2]<<'\n';
	j=sprintf(counterString,"%04d",finalFrame);		
	strcpy(fileName,argv[1]);
	strcat(fileName,counterString);
	strcat(fileName,".vtk");
	cout << fileName << ".\n";

// Set up colorbar
	vtkLookupTable *vortScale;
	vortScale = vtkLookupTable::New();
	vortScale->SetScaleToLinear();
	double range = max-min;
	int numColors = 100;
	vortScale->SetRange(min,max);
	vortScale->SetNumberOfTableValues(numColors);
	double colors[3][3] = { 0.0, 0.0, 1.0, //blue
							1.0, 1.0, 1.0, //white
							1.0, 0.0, 0.0 };//red
	double r,g,b;
	for (int i=0;i<numColors;i++){
		double val = min + ((double)i/numColors)*range;
		for (int j=0;j<2;j++){
			double currFloor = min + ((double)j/2.0) * range;
			double currCeil = min + ( (double)(j+1)/2.0)*range;
			if ( (val >= currFloor) && (val <=currCeil)){
				double currFraction = ( val - currFloor)/(currCeil-currFloor);
				r = colors[j][0]*(1.0-currFraction) + colors[j+1][0]*currFraction;
				g = colors[j][1]*(1.0-currFraction) + colors[j+1][1]*currFraction;
				b = colors[j][2]*(1.0-currFraction) + colors[j+1][2]*currFraction;
			}
		}
		vortScale->SetTableValue(i,r,g,b);
	}
	vortScale->Build();
	
// Reader loop	
	vtkPolyDataReader *sphereData[finalFrame];
	for ( int frameJ=0; frameJ<=finalFrame; frameJ++){
		
		sphereData[frameJ] = vtkPolyDataReader::New();
		
		j=sprintf(counterString,"%04d",frameJ);		
		strcpy(fileName,argv[1]);
		strcat(fileName,counterString);
		strcat(fileName,".vtk");
		cout << " Loading:"<< fileName << "...\n";
		
		sphereData[frameJ]->SetFileName(fileName);
		sphereData[frameJ]->SetScalarsName(scalarsDataName);
		
		}
	// Determine number of panels from data	
	int nActive = sphereData[0]->GetOutput()->GetNumberOfCells(); 
	// Write plot annotation string 
	char textString[64];
	j=sprintf(textString,"t = %7.4f, N = %d",0.000*dt,nActive);

	// Set up VTK Sphere plot pipeline
	vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
		sphereMapper->SetScalarRange(min,max);
		sphereMapper->SetLookupTable(vortScale);
	vtkActor *sphereActor = vtkActor::New();
		sphereActor->SetMapper(sphereMapper);
		if (edgesOn){sphereActor->GetProperty()->EdgeVisibilityOn();
			sphereActor->GetProperty()->SetEdgeColor(0.0,0.0,0.0);}
	// Set up VTK 3D axes pipeline	
// 	vtkAxesActor *axes = vtkAxesActor::New();
// 		axes->SetShaftTypeToLine();
// 		axes->SetXAxisLabelText("x");
// 		axes->SetYAxisLabelText("y");
// 		axes->SetZAxisLabelText("z");
// 		axes->SetTotalLength(1.5,1.5,1.5);
// 		axes->GetXAxisShaftProperty()->SetColor(0.6,0.6,0.6);
// 		axes->GetYAxisShaftProperty()->SetColor(0.6,0.6,0.6);
// 		axes->GetZAxisShaftProperty()->SetColor(0.6,0.6,0.6);
// 		axes->GetXAxisTipProperty()->SetColor(0.6,0.6,0.6);
// 		axes->GetYAxisTipProperty()->SetColor(0.6,0.6,0.6);
// 		axes->GetZAxisTipProperty()->SetColor(0.6,0.6,0.6);
// 	vtkCaptionActor2D *xAxisLabel = axes->GetXAxisCaptionActor2D();
// 		xAxisLabel->GetProperty()->SetColor(0.8,0.8,0.8);
// 	vtkCaptionActor2D *yAxisLabel = axes->GetYAxisCaptionActor2D();
// 		yAxisLabel->GetProperty()->SetColor(0.8,0.8,0.8);
// 	vtkCaptionActor2D *zAxisLabel = axes->GetZAxisCaptionActor2D();
// 		zAxisLabel->GetProperty()->SetColor(0.8,0.8,0.8);
	// Set up VTK scalar bar pipeline
	vtkScalarBarActor *scalarBar = vtkScalarBarActor::New();
		scalarBar->SetTitle(scalarName);
		scalarBar->SetOrientationToVertical();
		scalarBar->SetPosition(0.025,0.1);
		scalarBar->SetPosition2(0.9,0.9);
		scalarBar->SetWidth(0.2);
		scalarBar->SetLookupTable(vortScale);
	vtkTextProperty *scalarTextProperty=scalarBar->GetLabelTextProperty();
		//scalarTextProperty->SetColor(0.8,0.8,0.8);
		scalarTextProperty->SetColor(0,0,0);
	vtkTextProperty *scalarTitleTextProperty =scalarBar->GetTitleTextProperty();
		//scalarTitleTextProperty->SetColor(0.8,0.8,0.8);
		scalarTitleTextProperty->SetColor(0,0,0);
	
	// Set up VTK text annotation pipeline
	vtkTextActor *textActor = vtkTextActor::New();
		textActor->SetTextScaleModeToProp();
		textActor->SetPosition(0.1,0.1);
		textActor->SetPosition2(0.8,0.3);
		textActor->SetInput(textString);
	vtkTextProperty *annotationTextProperty=textActor->GetTextProperty();
		//annotationTextProperty->SetColor(0.8,0.8,0.8);
		annotationTextProperty->SetColor(0,0,0);
	// Set up camera
	vtkCamera *camera = vtkCamera::New();
		camera->SetPosition(5,0,2);
		camera->SetFocalPoint(0,0,0);
		camera->SetViewUp(0,0,1);
	
	// Set up renderer
	vtkRenderer *renderer = vtkRenderer::New();
	//	renderer->AddActor(axes);
		renderer->AddActor(sphereActor);
		renderer->AddActor(scalarBar);
		renderer->AddActor(textActor);
		renderer->SetActiveCamera(camera);
		renderer->SetBackground(0.9,0.9,0.9);
		//renderer->SetBackground(0,0,0);
	vtkRenderWindow *renWin = vtkRenderWindow::New();
		renWin->AddRenderer(renderer);
		renWin->SetSize(600,600);
	
	// Set up 3D to 2D filter and JPG writer
	vtkWindowToImageFilter *win2im = vtkWindowToImageFilter::New();
		win2im->SetInput(renWin);
	vtkJPEGWriter *writer = vtkJPEGWriter::New();
	cout << "...Loading complete.\n";
	cout << "Rendering...\n";
	// Render loop	
	for ( int frameJ=0; frameJ<=finalFrame; frameJ++){
		sphereMapper->SetInput(sphereData[frameJ]->GetOutput());
		sphereMapper->Update();
		scalarBar->SetLookupTable(sphereMapper->GetLookupTable());
		nActive = sphereData[frameJ]->GetOutput()->GetNumberOfCells();
		j=sprintf(textString,"t = %7.4f, N = %d",frameJ*dt,nActive);
		textActor->SetInput(textString);
		// Prepare output file
		j=sprintf(counterString,"%04d",frameJ);		
		strcpy(jpgFileName,jpgFileRoot);
		strcat(jpgFileName,counterString);
		strcat(jpgFileName,".jpg");
		// Render scene
		renWin->Render();
		// Write to JPG file
		win2im->Modified();
		writer->SetInput(win2im->GetOutput());
		writer->SetFileName(jpgFileName);
		writer->Write();
		
		if ( frameJ == finalFrame/4){cout << "... 25\% done.\n";}
		if ( frameJ == finalFrame/2){cout << "... 50\% done.\n";}
		if ( frameJ == 3*finalFrame/4){cout << "... 75\% done.\n";}
		}
	cout << "Rendering complete.\n";
	// Clean up
	for ( int frameJ=0;frameJ<=finalFrame;frameJ++){
		sphereData[frameJ]->Delete();
		}
	win2im->Delete();
	writer->Delete();
	sphereActor->Delete();
	sphereMapper->Delete();
//	axes->Delete();
	scalarBar->Delete();
	camera->Delete();
	renderer->Delete();
	renWin->Delete();
	textActor->Delete();
	vortScale->Delete();
};
