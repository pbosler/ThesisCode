// Pete Bosler
// February, 2013
// Animate .vtk relative vorticity plot initial latitude
// filename : vtkPlot2Panel.cpp
// Last updated : 2-16-2013
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

void GetColorForValue_BlueWhiteRed(double val, double &r, double &g, double &b, double &min, double &max);
void GetColorForValue_ColdToHot(double val, double &r, double &g, double &b, double &min, double &max);
void GetColorForValue_JetLat(double val, double &r, double &g, double &b, double &min, double &max);

// GOALS : one panel showing the whole sphere with relative vorticity
//		   one panel showing the adaptive mesh

int main(int argc, const char *argv[]){

	
	/* CHANGE DATA PARAMETERS HERE */
	/*	TO DO : this should happen at run time */
		double dt = 0.0025;
		double relVortMin = -12.0;
		double relVortMax = 12.0;
		char scalarsDataName[16] = "relVortPanel"; // must match .vtk file
		char scalarsDataTitle[64] = "Relative Vorticity";
		char scalarsDataName2[16] = "Tracer1";
		char scalarsDataTitle2[64] = "Initial Latitude";
		char titleString1[36] = "Perturbed Jet";
		char titleString2[56];
	/* END CHANGE DATA PARAMETERS */
	
	static char vtkFileRoot[128];
	static int finalFrame=0;
//	if ( argc != 3){
// 		cout << "Enter file root for .vtk files : \n";
// 		cin >> vtkFileRoot;
// 		cout << "\nEnter number of frames : \n";
// 		cin >> finalFrame;
//	}
// 	else{
		strcpy(vtkFileRoot,argv[1]);
		istringstream ss(argv[2]);
		if ( !(ss>>finalFrame) ) cerr << "Invalid final frame number.";
// 	}
	
	char counterString[4];
	int j;
	j = sprintf(counterString,"%04d",0);
	char vtkFileName[256];
	strcpy(vtkFileName,vtkFileRoot);
	strcat(vtkFileName,counterString);
	strcat(vtkFileName,".vtk");
	
	cout << "Preparing plot data ... \n";
	cout << "... reading " << finalFrame << " files beginning with : \n";
	cout << vtkFileName << "\n";
	
	// set up vtk file readers
	vtkPolyDataReader *sphereData[finalFrame+1];
	vtkPolyDataReader *sphereData2[finalFrame+1];
	for (int frameJ=0; frameJ<=finalFrame; frameJ++){
		// Set up input file names
		j=sprintf(counterString,"%04d",frameJ);
		strcpy(vtkFileName,vtkFileRoot);
		strcat(vtkFileName,counterString);
		strcat(vtkFileName,".vtk");
		// initialize readers
		sphereData[frameJ]=vtkPolyDataReader::New();
		sphereData[frameJ]->SetFileName(vtkFileName);	
		sphereData[frameJ]->SetScalarsName(scalarsDataName);
		
		sphereData2[frameJ]=vtkPolyDataReader::New();
		sphereData2[frameJ]->SetFileName(vtkFileName);
		sphereData2[frameJ]->SetScalarsName(scalarsDataName);
	}
	cout << "... done.\n";
	
	// Set up colorbars
	cout << "... building RelVort colorbar ...\n";
	vtkLookupTable *vortScale=vtkLookupTable::New();
	vortScale->SetScaleToLinear();
	int numColors = 512;
	double range = relVortMax - relVortMin;
	vortScale->SetRange(relVortMin,relVortMax);
	vortScale->SetNumberOfTableValues(numColors);
	double r, g, b, val;
	for (int i=0;i<numColors;i++){
		val = relVortMin + ((double)i/numColors)*range;
		GetColorForValue_BlueWhiteRed(val,r,g,b,relVortMin,relVortMax);
		vortScale->SetTableValue(i,r,g,b);
	}
	vortScale->Build();
	cout << "... done.\n";
	
	cout << "... building Flow Map colorbar ...\n";
	vtkLookupTable *latScale=vtkLookupTable::New();
	latScale->SetScaleToLinear();
	double latMin = -1.5708;
	double latMax = 1.5708;
	range = latMax-latMin;
	latScale->SetRange(latMin,latMax);
	latScale->SetNumberOfTableValues(numColors);
	for (int i=0;i<numColors;i++){
		val = latMin + ((double)i/numColors)*range;
		GetColorForValue_JetLat(val,r,g,b,latMin,latMax);
		latScale->SetTableValue(i,r,g,b);
	}
	latScale->Build();
	cout << "... done.\n";
	
	// Set up vtk pipeline 
	cout << "... building VTK graphics pipeline... \n";
	vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
		sphereMapper->SetScalarRange(relVortMin,relVortMax);
		sphereMapper->SetLookupTable(vortScale);
		//sphereMapper->SetScalarRange(latMin,latMax);
		//sphereMapper->SetLookupTable(latScale);
	vtkPolyDataMapper *sphereMapper2 = vtkPolyDataMapper::New();
		//sphereMapper2->SetScalarRange(latMin,latMax);
		//sphereMapper2->SetLookupTable(latScale);
		sphereMapper2->SetScalarRange(relVortMin,relVortMax);
		sphereMapper2->SetLookupTable(vortScale);
	vtkActor *sphereActor = vtkActor::New();
		sphereActor->SetMapper(sphereMapper);
	vtkActor *sphereActor2 = vtkActor::New();
		sphereActor2->SetMapper(sphereMapper2);
		sphereActor2->GetProperty()->EdgeVisibilityOn();
		sphereActor2->GetProperty()->SetEdgeColor(0.0,0.0,0.0);
	vtkScalarBarActor *colorbar = vtkScalarBarActor::New();
		colorbar->SetTitle(scalarsDataTitle);
		colorbar->GetTitleTextProperty()->SetColor(0.0,0.0,0.0);
		colorbar->SetOrientationToVertical();
		colorbar->SetPosition(0.025,0.1);
		colorbar->SetPosition2(0.9,0.9);
		colorbar->SetWidth(0.2);
		colorbar->SetLookupTable(vortScale);
		colorbar->GetLabelTextProperty()->SetColor(0.0,0.0,0.0);
	vtkScalarBarActor *colorbar2=vtkScalarBarActor::New();
		colorbar2->SetTitle(scalarsDataTitle2);
		colorbar2->SetOrientationToVertical();
		colorbar2->SetPosition(0.025,0.1);
		colorbar2->SetPosition2(0.9,0.9);
		colorbar2->SetWidth(0.2);
		colorbar2->SetLookupTable(latScale);
		colorbar2->GetLabelTextProperty()->SetColor(1.0,1.0,1.0);
	vtkTextActor *textActor1 = vtkTextActor::New();
		textActor1->SetTextScaleModeToProp();
		textActor1->SetPosition(0.1,0.1);
		textActor1->SetPosition2(0.8,0.3);
		textActor1->GetTextProperty()->SetColor(0.0,0.0,0.0);
		textActor1->SetInput(titleString1);
	vtkTextActor *textActor2 = vtkTextActor::New();
		textActor2->SetTextScaleModeToProp();
		textActor2->SetPosition(0.1,0.1);
		textActor2->SetPosition2(0.8,0.3);
		textActor2->GetTextProperty()->SetColor(0.0,0.0,0.0);
	vtkCamera *camera1 = vtkCamera::New();
		camera1->SetPosition(2.88,-1.44,2.64);
		camera1->SetFocalPoint(0.0,0.0,0.0);
		camera1->SetViewUp(-0.5,0.3,0.75);
	vtkCamera *camera2 = vtkCamera::New();
		camera2->SetPosition(2.1,-1.0,1.9);
		camera2->SetFocalPoint(0.0,0.0,0.0);
		camera2->SetViewUp(-0.5,0.3,0.75);
	vtkRenderer *renderer1 = vtkRenderer::New();
		renderer1->AddActor(sphereActor);
		renderer1->AddActor(colorbar);
		renderer1->AddActor(textActor2);
		renderer1->SetActiveCamera(camera1);
		renderer1->SetBackground(0.9,0.9,0.9);
		renderer1->SetViewport(0.0,0.0,0.5,1.0);
	vtkRenderer *renderer2 = vtkRenderer::New();
		renderer2->AddActor(sphereActor2);
		//renderer2->AddActor(colorbar2);
		renderer2->AddActor(textActor1);
		renderer2->SetActiveCamera(camera2);
		renderer2->SetBackground(0.9,0.9,0.9);
		renderer2->SetViewport(0.5,0.0,1.0,1.0);
	vtkRenderWindow *renWin = vtkRenderWindow::New();
		renWin->AddRenderer(renderer1);
		renWin->AddRenderer(renderer2);
		renWin->SetSize(1200,600);
	vtkWindowToImageFilter *win2im = vtkWindowToImageFilter::New();
		win2im->SetInput(renWin);
	vtkJPEGWriter *writer = vtkJPEGWriter::New();
		char jpgFileRoot[256];
		char jpgFileName[256];
		strcpy(jpgFileRoot,"jpgOut/jetAMR2panel_");
	cout << "... done.\n";
	
	// RENDERING LOOP
	cout << "... RENDERING ...\n";
	int nActive = 0;
	for (int frameJ = 0; frameJ<=finalFrame; frameJ++){
		sphereMapper->SetInput(sphereData[frameJ]->GetOutput());
		sphereMapper2->SetInput(sphereData2[frameJ]->GetOutput());
		
		j = sprintf(titleString2,"t = %7.4f, N = %6d",frameJ*dt,nActive);
		textActor2->SetInput(titleString2);
		
		j = sprintf(counterString,"%04d",frameJ);
		strcpy(jpgFileName,jpgFileRoot);
		strcat(jpgFileName,counterString);
		strcat(jpgFileName,".jpg");
		
		renWin->Render();
		
		win2im->Modified();
		writer->SetInput(win2im->GetOutput());
		writer->SetFileName(jpgFileName);
		writer->Write();
		
		nActive = sphereData[frameJ]->GetOutput()->GetNumberOfCells();
		
		if ( frameJ == finalFrame/4){cout << "... 25\% done.\n";}
		if ( frameJ == finalFrame/2){cout << "... 50\% done.\n";}
		if ( frameJ == 3*finalFrame/4){cout << "... 75\% done.\n";}
	}	
	// Plot last frame
	sphereMapper->SetInput(sphereData[finalFrame]->GetOutput());
	j = sprintf(counterString,"%04d",finalFrame+1);
	strcpy(jpgFileName,jpgFileRoot);
	strcat(jpgFileName,counterString);
	strcat(jpgFileName,".jpg");
	renWin->Render();
	win2im->Modified();
	writer->SetInput(win2im->GetOutput());
	writer->SetFileName(jpgFileName);
	writer->Write();
	
	// Clean up
	win2im->Delete();
	writer->Delete();
	renWin->Delete();
	renderer2->Delete();
	renderer1->Delete();
	camera2->Delete();
	camera1->Delete();
	textActor1->Delete();
	textActor2->Delete();
	colorbar->Delete();
	colorbar2->Delete();
	sphereActor->Delete();
	sphereActor2->Delete();
	sphereMapper->Delete();
	sphereMapper2->Delete();
	vortScale->Delete();
	latScale->Delete();
	for (int frameJ=0; frameJ<=finalFrame; frameJ++){
		sphereData[frameJ]->Delete();
		sphereData2[frameJ]->Delete();
	}
};
void GetColorForValue_BlueWhiteRed(double val, double &r, double &g, double &b, double &min,double &max)
{
	const int numColors = 3;
	const double range = max - min;
	const double colors[numColors][3] =  { 0.0, 0.0, 1.0, //blue
							1.0, 1.0, 1.0, //white
							1.0, 0.0, 0.0 };//red
	if ( val <= min )
	{
		r = colors[0][0];
		g = colors[0][1];
		b = colors[0][2];
	}
	else if ( val > max)
	{
		r = colors[numColors-1][0];
		g = colors[numColors-1][1];
		b = colors[numColors-1][2];
	}
	else
	{							
		for (int j=0;j<numColors;j++)
		{
			double floorj = min + ((double)j/(numColors-1))*range;
			double ceilj = min + ((double)(j+1)/(numColors-1))*range;
			if ( (val > floorj) && (val <= ceilj) ){
				double frac = ( val - floorj)/(ceilj-floorj);
				r = colors[j][0]*(1.0-frac) + colors[j+1][0]*frac;
				g = colors[j][1]*(1.0-frac) + colors[j+1][1]*frac;
				b = colors[j][2]*(1.0-frac) + colors[j+1][2]*frac;
			}
		}
	}
};

void GetColorForValue_JetLat(double val, double &r, double &g, double &b, double &min, double &max)
{
	const int numColors = 9;
	const double range = max-min;
	const double colors[numColors][3] = {0.0, 1.0, 1.0, //cyan
										 0.0, 0.0, 1.0,//blue
						  				 0.8549, 0.4392157, 0.839216,//orchid
						   				 0.51765, 0.4392157, 1.0,//light slate blue
						   				 0.28235, 0.8196078, 0.8,//medium turqoise
										 0.0, 1.0, 0.0,//green
						   				 1.0, 0.843137, 0.0,//gold
						   				 1.0, 0.54902, 0.0,//dark orange
						   				 1.0, 1.0, 0.0};//yellow
	double frac;
	double floor;
	double ceil;
	double points[9] = {-1.5708, -0.7854, 0.0, 0.6, 0.7, 0.7854, 0.9, 1.0, 1.5708};
	int color0;
	
	double c1, c2, c3;
	double p1, p2, p3;
	int color1, color2, color3;
	for (int j=0; j < numColors-1; j++)
	{
		if (( val > points[j]) && ( val<=points[j+1]))
		{
			floor = points[j];
			ceil = points[j+1];
			color0 = j;
		}
	}

	/* Linear color weighting */
	if ( val <= points[0])
		{
			r = colors[0][0];
		 	g = colors[0][1];
		 	b = colors[0][2];
		 }
	else if ( val > points[numColors-1])
		{
		   r = colors[numColors-1][0];
		   g = colors[numColors-1][1];
		   b = colors[numColors-1][2];
		}
	else
		{
		frac = (val-floor)/(ceil-floor);
		r = colors[color0][0]*(1.0-frac) + colors[color0+1][0]*frac;
		g = colors[color0][1]*(1.0-frac) + colors[color0+1][1]*frac;
		b = colors[color0][2]*(1.0-frac) + colors[color0+1][2]*frac;
		}

};

void GetColorForValue_ColdToHot(double val, double &r, double &g, double &b, double &min, double &max)
{
	const int numColors = 7;
	const double range = max - min;
	const double colors[numColors][3] = { 0.87843, 1.0, 1.0, //light cyan
							0.0, 1.0, 1.0, // cyan
							0.0, 0.0, 1.0, // blue
							0.0, 0.0, 0.0, // black
							1.0, 0.0, 0.0, // red
							1.0, 1.0, 0.0, // yellow
							1.0, 1.0, 0.87843}; // light yellow
	for (int j=0;j<numColors;j++){
		double floorj = min + ((double)j/(numColors-1))*range;
		double ceilj = min +((double)(j+1)/(numColors-1))*range;
		if ( (val > floorj) && (val <= ceilj) ) {
			double frac = ( val-floorj)/(ceilj-floorj);
			r = colors[j][0]*(1.0-frac) + colors[j+1][0]*frac;
			g = colors[j][1]*(1.0-frac) + colors[j+1][1]*frac;
			b = colors[j][2]*(1.0-frac) + colors[j+1][2]*frac;
		}
	}							

};