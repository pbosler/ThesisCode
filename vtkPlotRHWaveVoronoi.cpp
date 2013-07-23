// Pete Bosler
// April 17, 2012
// Animate .vtk relative vorticity plot initial latitude
// filename : vtkPlot2Panel.cpp
// Last updated : 9-1-2012
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
void GetColorForValue_MaizeBlue(double val, double &r, double &g, double &b, double &min, double &max);
//void GetColorForValue_LongitudeSectors(double val, double &r, double &g, double &b, double &min, double &max);

// GOALS : one panel showing the whole sphere with relative vorticity
//		   one panel showing a zoomed-in plot of the Voronoi panels and active particles with edges

int main(){
	//
	char voronoiFileRoot[128];
	char delaunayFileRoot[128];
	int finalFrame;
	cout << "Enter file root for Voronoi .vtk files : \n ";
	cin >> voronoiFileRoot;
	cout << "\n Enter file root for Delaunay .vtk files : \n";
	cin >> delaunayFileRoot;
	cout << "\n Enter number of frames :\n";
	cin >> finalFrame;
	
/* CHANGE DATA PARAMETERS HERE */
	double dt = 0.01;
	double relVortMin = -9.0;
	double relVortMax = 9.0;
	char vorScalarsDataName[16] = "relVortPanel"; // must match .vtk file
	char vorScalarsTitle[64] = "Relative Vorticity";
	char delScalarsDataName[16] = "relVortPanel"; // must match .vtk file
	char delScalarsTitle[64] = "Relative Vorticity";
/* END DATA PARAMETERS SECTION */
/* TO DO : READ DATA PARAMETERS AT RUN TIME */	

	char titleString1[64];
	char titleString2[36];
		
// Set up input filenames
	char counterString[4];
	int counter = 0;
	int j;
	char vorfileName[128];
	char delfileName[128];
	
	j = sprintf(counterString,"%04d",counter);
	strcpy(vorfileName,voronoiFileRoot);
	strcpy(delfileName,delaunayFileRoot);	
	strcat(vorfileName,counterString);
	strcat(delfileName,counterString);
	strcat(vorfileName,".vtk");
	strcat(delfileName,".vtk");
	cout << "Preparing plot data ... \n";
	cout << " ... reading "<< finalFrame << " files beginning with : \n" 
		<< vorfileName << " and " << delfileName << "\n";

// Set up VTK File Readers
	vtkPolyDataReader *DelaunayData[finalFrame+1];
	vtkPolyDataReader *VoronoiData[finalFrame+1];
	vtkPolyDataReader *DelaunayData2[finalFrame+1];
	vtkPolyDataReader *VoronoiData2[finalFrame+1];
	for (int frameJ=0; frameJ<=finalFrame; frameJ++){
		// initialize readers
		DelaunayData[frameJ]=vtkPolyDataReader::New();
		VoronoiData[frameJ]=vtkPolyDataReader::New();
		DelaunayData2[frameJ]=vtkPolyDataReader::New();
		VoronoiData2[frameJ]=vtkPolyDataReader::New();
		// set readers' input files
		j=sprintf(counterString,"%04d",frameJ);
		strcpy(vorfileName,voronoiFileRoot);
		strcat(vorfileName,counterString);
		strcat(vorfileName,".vtk");
		strcpy(delfileName,delaunayFileRoot);
		strcat(delfileName,counterString);
		strcat(delfileName,".vtk");
		
		DelaunayData[frameJ]->SetFileName(delfileName);
		VoronoiData[frameJ]->SetFileName(vorfileName);
		DelaunayData2[frameJ]->SetFileName(delfileName);
		VoronoiData2[frameJ]->SetFileName(vorfileName);
		
		// set readers' scalar data source
		//DelaunayData[frameJ]->SetScalarsName(delScalarsDataName);
		VoronoiData[frameJ]->SetScalarsName(vorScalarsDataName);
		//DelaunayData2[frameJ]->SetScalarsName("tracerPanel1");
		VoronoiData2[frameJ]->SetScalarsName("tracerPanel1");
	}	
	cout << "... done.\n";
	cout << "... building RelVort Colorbar ... \n";
// set up colorbars
	// Vorticity colorbar
	vtkLookupTable *vortScale=vtkLookupTable::New();
	vortScale->SetScaleToLinear();
	int numColors=512;
	double range=relVortMax - relVortMin;
	vortScale->SetRange(relVortMin,relVortMax);
	vortScale->SetNumberOfTableValues(numColors);
	double r, g, b;
	double val;
	for (int i=0; i<numColors; i++){
		val = relVortMin + ((double)i/numColors)*range;
		GetColorForValue_BlueWhiteRed(val,r,g,b,relVortMin,relVortMax);
		vortScale->SetTableValue(i,r,g,b);
	}
	vortScale->Build();
	
	// setup latitude colorbar
	vtkLookupTable *latScale=vtkLookupTable::New();
	double latMin = -1.5708;
	double latMax = 1.5708;
	range = latMax - latMin;
	latScale->SetScaleToLinear();
	latScale->SetRange(latMin,latMax);
	latScale->SetNumberOfTableValues(numColors);
	for (int i=0;i<numColors;i++){
		val = latMin + ((double)i/numColors)*range;
		GetColorForValue_MaizeBlue(val,r,g,b,latMin,latMax);
		latScale->SetTableValue(i,r,g,b);
	}
	latScale->Build();
	
	
	cout << "... done.\n";
	
// Setup vtk pipeline for left panel (whole sphere, relVort)
	cout << "... building vtk graphics pipeline...\n";
	vtkPolyDataMapper *delMapper1 = vtkPolyDataMapper::New();
	vtkPolyDataMapper *vorMapper1 = vtkPolyDataMapper::New();
		delMapper1->SetScalarRange(relVortMin,relVortMax);
		vorMapper1->SetScalarRange(relVortMin,relVortMax);
		//delMapper1->SetLookupTable(vortScale);
		vorMapper1->SetLookupTable(vortScale);
	vtkActor *delActor1=vtkActor::New();
	vtkActor *vorActor1=vtkActor::New();
		delActor1->SetMapper(delMapper1);
		delActor1->GetProperty()->SetRepresentationToPoints();
		delActor1->GetProperty()->SetPointSize(2.0);
		delActor1->GetProperty()->SetColor(0.0,0.0,0.0);
		vorActor1->SetMapper(vorMapper1);
	vtkScalarBarActor *relVortColorbar=vtkScalarBarActor::New();
		relVortColorbar->SetTitle(vorScalarsTitle);
		relVortColorbar->SetOrientationToVertical();
		relVortColorbar->SetPosition(0.025,0.1);
		relVortColorbar->SetPosition2(0.9,0.9);
		relVortColorbar->SetWidth(0.2);
		relVortColorbar->SetLookupTable(vortScale);
		vtkTextProperty *colorbarTextProperty1=relVortColorbar->GetLabelTextProperty();
			colorbarTextProperty1->SetColor(0,0,0);
	vtkTextActor *textActor1=vtkTextActor::New();
		textActor1->SetTextScaleModeToProp();
		textActor1->SetPosition(0.1,0.1);
		textActor1->SetPosition2(0.8,0.3);
		vtkTextProperty *annotationTextProperty=textActor1->GetTextProperty();
			annotationTextProperty->SetColor(0,0,0);
	vtkCamera *camera1=vtkCamera::New();
		camera1->SetPosition(5,0,2);
		camera1->SetFocalPoint(0,0,0);
		camera1->SetViewUp(0,0,1);
	

// Setup vtk pipeline for right panel (zoom in)	
	vtkPolyDataMapper *delMapper2=vtkPolyDataMapper::New();
	vtkPolyDataMapper *vorMapper2=vtkPolyDataMapper::New();
		vorMapper2->SetScalarRange(latMin,latMax);
		vorMapper2->SetLookupTable(latScale);
	vtkActor *delActor2=vtkActor::New();
		delActor2->GetProperty()->SetRepresentationToPoints();
		delActor2->GetProperty()->SetPointSize(4.0);
		delActor2->GetProperty()->SetColor(0.0,0.0,0.0);
	vtkActor *vorActor2=vtkActor::New();
		delActor2->SetMapper(delMapper2);
		vorActor2->SetMapper(vorMapper2);
		vorActor2->GetProperty()->EdgeVisibilityOn();
		vorActor2->GetProperty()->SetEdgeColor(0.0,0.0,0.0);
	vtkTextActor *textActor2=vtkTextActor::New();
		textActor2->SetTextScaleModeToProp();
		textActor2->SetPosition(0.1,0.1);
		textActor2->SetPosition(0.8,0.3);
		vtkTextProperty *annotationTextProperty2=textActor2->GetTextProperty();
			annotationTextProperty2->SetColor(0.0,0.0,0.0);	
	vtkScalarBarActor *colorbar2=vtkScalarBarActor::New();
		colorbar2->SetTitle("initial latitude");
		colorbar2->GetTitleTextProperty()->SetColor(0.0,0.0,0.0);
		colorbar2->GetLabelTextProperty()->SetColor(0.0,0.0,0.0);
		colorbar2->SetOrientationToVertical();
		colorbar2->SetPosition(0.025,0.1);
		colorbar2->SetPosition2(0.9,0.9);
		colorbar2->SetWidth(0.2);
		colorbar2->SetLookupTable(latScale);		
	

	cout << "... done. \n";
	

// Setup viewing window	
	cout << "... building vtk output ...\n";
	
	vtkRenderer *renderer1=vtkRenderer::New();
		renderer1->AddActor(delActor1);
		renderer1->AddActor(vorActor1);
		renderer1->AddActor(relVortColorbar);
		renderer1->AddActor(textActor2);
		renderer1->SetActiveCamera(camera1);
		renderer1->SetBackground(0.9,0.9,0.9);
		renderer1->SetViewport(0.0,0.0,0.5,1.0);
	
	vtkRenderer *renderer2=vtkRenderer::New();
		renderer2->AddActor(delActor2);
		renderer2->AddActor(vorActor2);
		renderer2->AddActor(colorbar2);
		renderer2->AddActor(textActor1);
		renderer2->SetActiveCamera(camera1);
		renderer2->SetBackground(0.9,0.9,0.9);
		renderer2->SetViewport(0.5,0.0,1.0,1.0);		
	
	vtkRenderWindow *renWin=vtkRenderWindow::New();
		renWin->AddRenderer(renderer1);
		renWin->AddRenderer(renderer2);
		renWin->SetSize(1200,600);
	vtkWindowToImageFilter *win2im=vtkWindowToImageFilter::New();
		win2im->SetInput(renWin);
	



// Set up output filenames and writer
	char jpgFileRoot[256];
	char jpgFileName[256];
	strcpy(jpgFileName,"jpgOut/rh4Voronoi_");
	strcpy(jpgFileRoot,jpgFileName);
	strcat(jpgFileName,counterString);
	strcat(jpgFileName,".jpg");	

	vtkJPEGWriter *writer=vtkJPEGWriter::New();
	
	cout << "... done.\n";
	
// Rendering loop
	cout << "Rendering ...\n";
	int nActive = 0;
	textActor1->SetInput(titleString1);
	for (int frameJ=0; frameJ<=finalFrame; frameJ++){	
		delMapper1->SetInput(DelaunayData[frameJ]->GetOutput());
		delMapper2->SetInput(DelaunayData2[frameJ]->GetOutput());
		vorMapper1->SetInput(VoronoiData[frameJ]->GetOutput());
		vorMapper2->SetInput(VoronoiData2[frameJ]->GetOutput());
		//relVortColorbar->SetLookupTable(vorMapper1->GetLookupTable);
		
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
		
		nActive=VoronoiData[frameJ]->GetOutput()->GetNumberOfCells();
		
		if ( frameJ == finalFrame/4){cout << "... 25\% done.\n";}
		if ( frameJ == finalFrame/2){cout << "... 50\% done.\n";}
		if ( frameJ == 3*finalFrame/4){cout << "... 75\% done.\n";}
	}
	// Plot last frame
	delMapper1->SetInput(DelaunayData[finalFrame]->GetOutput());
	delMapper2->SetInput(DelaunayData2[finalFrame]->GetOutput());
	vorMapper1->SetInput(VoronoiData[finalFrame]->GetOutput());
	vorMapper2->SetInput(VoronoiData2[finalFrame]->GetOutput());
	j = sprintf(titleString2,"t = %7.4f, N = %6d",finalFrame*dt,nActive);
	textActor2->SetInput(titleString2);
	j = sprintf(counterString,"%04d",finalFrame+1);
	strcpy(jpgFileName,jpgFileRoot);
	strcat(jpgFileName,counterString);
	strcat(jpgFileName,".jpg");
	renWin->Render();
	win2im->Modified();
	writer->SetInput(win2im->GetOutput());
	writer->SetFileName(jpgFileName);
	writer->Write();
	cout << "Rendering complete.\n";
	
// Cleanup
	for (int frameJ=0;frameJ<=finalFrame;frameJ++){
		DelaunayData[frameJ]->Delete();
		VoronoiData[frameJ]->Delete();
		DelaunayData2[frameJ]->Delete();
		VoronoiData2[frameJ]->Delete();
	}	
	win2im->Delete();
	writer->Delete();
	renWin->Delete();
	renderer2->Delete();
	renderer1->Delete();
	camera1->Delete();
	vorActor2->Delete();
	vorActor1->Delete();
	delActor2->Delete();
	delActor1->Delete();
	textActor2->Delete();
	textActor1->Delete();
	vorMapper1->Delete();
	vorMapper2->Delete();
	delMapper1->Delete();
	delMapper2->Delete();
	relVortColorbar->Delete();
	vortScale->Delete();
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

void GetColorForValue_MaizeBlue(double val, double &r, double &g, double &b, double &min,double &max)
{
	const int numColors = 3;
	const double range = max - min;
	const double colors[numColors][3] = { 0.0, 0.0, 0.4, // midnight blue
										  1.0, 1.0, 1.0, // white
										  1.0, 0.8, 0.2}; // yellow
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
// void GetColorForValue_LongitudeSectors(double val, double &r, double &g, double &b, double &min, double &max);
// {
// 	const int numcolors = 3;
// 	
// };
