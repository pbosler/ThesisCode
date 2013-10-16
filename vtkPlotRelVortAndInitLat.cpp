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

int main( int argc, const char *argv[]){
	// Check for 3 arguments (script name, plus 2 user arguments)
	if ( argc != 3){
		cout << "Usage : "<< argv[0] << " <jobtitle | frame total>\n";
		return 1;
	};
	char filename[128];
	char counterString[4];
	int j;
	int counter = 0;
	int finalFrame;
	char jpgFilename[256];
	char jpgFileRoot[256];
	char textString1[64];
	char textString2[64];
		
	/* USERS ALTER THIS SECTION */
	float dt = 0.01;	// time step size
	double min = -4.0; // minimum scalar value to be plotted // TO DO: Determine this from the vtk data file
	double max = 12.0;  // maximum scalar value to be plotted // TO DO: Determine this from the vtk data file
	char scalarsDataName1[16] = "relVortPanel"; // vtk data name (must match .vtk file)
	char scalarsDataName2[16] = "Tracer1" ;// vtk data name (must match .vtk file)
	char scalarTitle1[64] = "Relative Vorticity"; // scalar label for plot 1
	char scalarTitle2[64] = "Initial latitude"; // scalar label for plot 2
	char titleString[64] = " "; // title for both plots
	bool panel2LogScale = false; // use logarithmic color scale 
	bool edgesOn = true;
	/* END USERS ALTERATION SECTION */
	
	int nActiveArray[21] = { 3900, 3936, 3984, 4056, 4080,
							 4092, 4656, 5832, 7008, 7968,
							 8760, 9432, 10666, 10128, 10452,
							 10680, 11076, 11364, 11556, 11736,
							 11892};
	int remeshCounter = 1;
	int remeshNumber = 0;
	int remeshInterval = 20;
	int nActive = 0;
	cout << "initial nActive = " << nActive << "\n";
	
	
	// setup input filenames
	j = sprintf(counterString,"%04d",counter);
	strcpy(filename,argv[1]);
	strcat(filename,counterString);
	strcat(filename,".vtk");
	// setup output filenames
	strcpy(jpgFilename,"jpgOut/twoPanelMovie_");
	strcpy(jpgFileRoot,jpgFilename);
	strcat(jpgFilename,counterString);
	strcat(jpgFilename,".jpg");
	
	istringstream ss(argv[2]);
	if ( !(ss>>finalFrame) ) cerr << "Invalid final frame number "<<argv[2]<<'\n';
	cout << "Preparing plot data ... \n";
	cout << "... reading "<< finalFrame << " files beginning with "<< filename << '\n';
	
	// File reader arrays
	vtkPolyDataReader *sphereData[finalFrame+1];
	vtkPolyDataReader *sphereData2[finalFrame+1];
	for ( int frameJ=0; frameJ<=finalFrame; frameJ++){
		sphereData[frameJ]=vtkPolyDataReader::New();
		sphereData2[frameJ]=vtkPolyDataReader::New();
		j=sprintf(counterString,"%04d",frameJ);
		strcpy(filename,argv[1]);
		strcat(filename,counterString);
		strcat(filename,".vtk");
		
		sphereData[frameJ]->SetFileName(filename);
		sphereData[frameJ]->SetScalarsName(scalarsDataName1);
		
		sphereData2[frameJ]->SetFileName(filename);
		sphereData2[frameJ]->SetScalarsName(scalarsDataName2);
	}
	// Determine number of panels from data
	//nActive = nActiveArray[0];
	//int nActive=sphereData[0]->GetOutput()->GetNumberOfCells();
	j = sprintf(textString1,"t = %7.4f, N = %d",0.0*dt,nActive);
	
	// set up vorticity colorbar
	vtkLookupTable *vortScale=vtkLookupTable::New();
	vortScale->SetScaleToLinear();
	int numColors = 512;
	double range = max-min;
	vortScale->SetRange(min,max);
	vortScale->SetNumberOfTableValues(numColors);
	double r, g, b;
	for (int i=0;i<numColors;i++){
		double val = min + ((double)i/numColors)*range;
		GetColorForValue_BlueWhiteRed(val,r,g, b,min,max);
		vortScale->SetTableValue(i,r,g,b);
	}
	vortScale->Build();
	
	// setup latitude colorbar
	vtkLookupTable *latScale=vtkLookupTable::New();
	int numColors2 = 512;
	double min2 = -1.5708;
	double max2 = 1.5708;
	latScale->SetScaleToLinear();
	double range2 = max2-min2;
	latScale->SetRange(min2,max2);
	latScale->SetNumberOfTableValues(numColors2);
	for (int i=0;i<numColors2;i++){
		double val2 = min2 + ((double)i/numColors2)*range2;
		GetColorForValue_JetLat(val2,r,g,b,min2,max2);
		latScale->SetTableValue(i,r,g,b);
	}
	latScale->Build();
	
	// Setup vtk pipeline for left panel
	vtkPolyDataMapper *sphereMapper1 = vtkPolyDataMapper::New();
		sphereMapper1->SetScalarRange(min,max);
		sphereMapper1->SetLookupTable(vortScale);
	vtkActor *sphereActor1 = vtkActor::New();
		sphereActor1->SetMapper(sphereMapper1);
		//if (edgesOn){sphereActor1->GetProperty()->EdgeVisibilityOn();}
	vtkScalarBarActor *scalarBar1 = vtkScalarBarActor::New();
		scalarBar1->SetTitle(scalarTitle1);
		scalarBar1->SetOrientationToVertical();
		scalarBar1->SetPosition(0.025,0.1);
		scalarBar1->SetPosition2(0.9,0.9);
		scalarBar1->SetWidth(0.2);
		scalarBar1->SetLookupTable(vortScale);
		vtkTextProperty *scalarTextProperty1=scalarBar1->GetLabelTextProperty();
			//scalarTextProperty1->SetColor(0.8,0.8,0.8);
			scalarTextProperty1->SetColor(0,0,0);
		vtkTextProperty *scalarTitleTextProperty1 =scalarBar1->GetTitleTextProperty();
			//scalarTitleTextProperty1->SetColor(0.8,0.8,0.8);
			scalarTitleTextProperty1->SetColor(0,0,0);
	vtkTextActor *textActor1=vtkTextActor::New();
		textActor1->SetTextScaleModeToProp();
		textActor1->SetPosition(0.1,0.1);
		textActor1->SetPosition2(0.8,0.3);
		vtkTextProperty *annotationTextProperty1=textActor1->GetTextProperty();
		//annotationTextProperty1->SetColor(0.8,0.8,0.8);
		annotationTextProperty1->SetColor(0,0,0);
	
	
	// Setup vtk pipeline for right panel
	vtkPolyDataMapper *sphereMapper2 = vtkPolyDataMapper::New();
		sphereMapper2->SetLookupTable(latScale);
		sphereMapper2->SetScalarRange(min2,max2);
	vtkActor *sphereActor2 = vtkActor::New();
		sphereActor2->SetMapper(sphereMapper2);
		if (edgesOn){
			sphereActor2->GetProperty()->EdgeVisibilityOn();
			//sphereActor2->GetProperty()->SetEdgeColor(0.098039,0.098039,0.4392156);
			sphereActor2->GetProperty()->SetEdgeColor(0.0, 1.0, 0.0);
			}
//TO DO: Set log scale	if ( panel2LogScale ) sphereMapper2->GetLookupTable()->SetScaleToLog10();
	vtkScalarBarActor *scalarBar2 = vtkScalarBarActor::New();
		scalarBar2->SetTitle(scalarTitle2);
		scalarBar2->SetOrientationToVertical();
		scalarBar2->SetPosition(0.025,0.1);
		scalarBar2->SetPosition2(0.9,0.9);
		scalarBar2->SetWidth(0.2);
		scalarBar2->SetLookupTable(latScale);
		vtkTextProperty *scalarTextProperty2=scalarBar2->GetLabelTextProperty();
			//scalarTextProperty2->SetColor(0.8,0.8,0.8);
			scalarTextProperty2->SetColor(0,0,0);
		vtkTextProperty *scalarTitleTextProperty2 =scalarBar2->GetTitleTextProperty();
			//scalarTitleTextProperty2->SetColor(0.8,0.8,0.8);
			scalarTitleTextProperty2->SetColor(0,0,0);
	vtkTextActor *textActor2=vtkTextActor::New();
		textActor2->SetTextScaleModeToProp();
		textActor2->SetPosition(0.1,0.1);
		textActor2->SetPosition2(0.8,0.3);
		vtkTextProperty *annotationTextProperty2=textActor2->GetTextProperty();
		//annotationTextProperty2->SetColor(0.8,0.8,0.8);
		annotationTextProperty2->SetColor(0,0,0);
	
	vtkTextActor *titleTextActor = vtkTextActor::New();		
		titleTextActor->SetTextScaleModeToProp();
		//titleTextActor->SetPosition(0.1,0.1);
		titleTextActor->SetPosition2(1.0,5.0);
		vtkTextProperty *annotationTextProperty3=titleTextActor->GetTextProperty();
		//annotationTextProperty3->SetColor(0.8,0.8,0.8);
		annotationTextProperty3->SetColor(0,0,0);
		titleTextActor->SetInput(titleString);
		
	// setup camera for both panels
	// Set up camera
	vtkCamera *camera = vtkCamera::New();
		camera->SetPosition(2.6,-3.9,2.6);
		//camera->SetPosition(5,0,2);
		camera->SetFocalPoint(0,0,0);
		camera->SetViewUp(0,0,1);
		
	
	// setup renderers for both panels
	vtkRenderer *renderer1 = vtkRenderer::New();
		renderer1->AddActor(sphereActor1);
		renderer1->AddActor(scalarBar1);
		renderer1->AddActor(textActor1);
		//renderer1->AddActor(titleTextActor);
		renderer1->SetActiveCamera(camera);
		//renderer1->SetBackground(0,0,0);
		renderer1->SetBackground(0.9,0.9,0.9);
		renderer1->SetViewport(0.0,0.0,0.5,1.0);
	vtkRenderer *renderer2 = vtkRenderer::New();
		renderer2->AddActor(sphereActor2);
		renderer2->AddActor(scalarBar2);
		renderer2->AddActor(textActor2);
		renderer2->SetActiveCamera(camera);
		//renderer2->SetBackground(0,0,0);		
		renderer2->SetBackground(0.9,0.9,0.9);		
		renderer2->SetViewport(0.5,0.0,1.0,1.0);
	vtkRenderWindow *renWin=vtkRenderWindow::New();
		renWin->AddRenderer(renderer1);
		renWin->AddRenderer(renderer2);
		renWin->SetSize(1200,600);
		
	// Setup 3D to 2D filter and JPEG writer
	vtkWindowToImageFilter *win2im = vtkWindowToImageFilter::New();
		win2im->SetInput(renWin);
	vtkJPEGWriter *writer = vtkJPEGWriter::New();
	
	// loop over all file readers
	for ( int frameJ = 0; frameJ <=finalFrame; frameJ++){
		sphereMapper1->SetInput(sphereData[frameJ]->GetOutput());
		sphereMapper2->SetInput(sphereData2[frameJ]->GetOutput());
		scalarBar1->SetLookupTable(sphereMapper1->GetLookupTable());
		scalarBar2->SetLookupTable(sphereMapper2->GetLookupTable());
		
			
		j = sprintf(textString1,"t = %7.4f, N = %6d",frameJ*dt,nActive);
		j = sprintf(textString2,"Perturbed Jet (no rotation)");
		textActor1->SetInput(textString1);
		textActor2->SetInput(textString2);
		
		j=sprintf(counterString,"%04d",frameJ);
		strcpy(jpgFilename,jpgFileRoot);
		strcat(jpgFilename,counterString);
		strcat(jpgFilename,".jpg");
		
		renWin->Render();
		
		win2im->Modified();
		writer->SetInput(win2im->GetOutput());
		writer->SetFileName(jpgFilename);
		writer->Write();
		
		nActive=sphereData2[frameJ]->GetOutput()->GetNumberOfCells();
		
		remeshCounter++;
		if ( remeshCounter%remeshInterval == 0 ) {
			remeshNumber++;
			cout << "Remesh" << remeshNumber<<" new nActive should be " << nActiveArray[remeshNumber] <<"\n";
			cout << " nActive is " << nActive <<"\n";
			}
//		nActive = nActiveArray[remeshNumber];
//		nActive = 20480;
		
		if ( frameJ == finalFrame/4){cout << "... 25\% done.\n";}
		if ( frameJ == finalFrame/2){cout << "... 50\% done.\n";}
		if ( frameJ == 3*finalFrame/4){cout << "... 75\% done.\n";}
	}
	// Plot Last Frame
	sphereMapper1->SetInput(sphereData[finalFrame]->GetOutput());
	sphereMapper2->SetInput(sphereData2[finalFrame]->GetOutput());
	scalarBar1->SetLookupTable(sphereMapper1->GetLookupTable());
	scalarBar2->SetLookupTable(sphereMapper2->GetLookupTable());
	
		
	j = sprintf(textString1,"t = %7.4f, N = %6d",finalFrame*dt,nActive);
	j = sprintf(textString2,"Steady Jet (no rotation)");
	textActor1->SetInput(textString1);
	textActor2->SetInput(textString2);
	
	j=sprintf(counterString,"%04d",finalFrame+1);
	strcpy(jpgFilename,jpgFileRoot);
	strcat(jpgFilename,counterString);
	strcat(jpgFilename,".jpg");
	
	renWin->Render();
	
	win2im->Modified();
	writer->SetInput(win2im->GetOutput());
	writer->SetFileName(jpgFilename);
	writer->Write();
	
	cout << "RemeshNumber = "<< remeshNumber << "\n";	
	cout << "Rendering complete.\n";
	
	// Cleanup
	for (int frameJ = 0; frameJ<=finalFrame;frameJ++){
		sphereData[frameJ]->Delete();
		sphereData2[frameJ]->Delete();
	}
	win2im->Delete();
	writer->Delete();
	renWin->Delete();
	renderer2->Delete();
	renderer1->Delete();
	camera->Delete();
	scalarBar2->Delete();
	scalarBar1->Delete();
	textActor2->Delete();
	textActor1->Delete();
	sphereMapper2->Delete();
	sphereMapper1->Delete();
	sphereActor2->Delete();
	sphereActor1->Delete();
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
	
	
// 	if ( val <= points[1] )
// 	{
// 		color1 = 0;
// 		color2 = 1;
// 		color3 = 2;
// 		floor = points[0];
// 		ceil = points[1];
// 		p1 = points[0];
// 		p2 = points[1];
// 		p3 = points[2];
// 	}
// 	else if ( val <= points[numColors-2] )
// 		for (int jj = 1;jj<=numColors-2;jj++)
// 		{
// 			if ( (val > points[jj]) && ( val <= points[jj+1]))
// 			{
// 				color1 = jj-1;
// 				color2 = jj;
// 				color3 = jj+1;				
// 				p1 = points[jj-1];
// 				p2 = points[jj];
// 				p3 = points[jj+1];
// 				
// 				floor = points[jj];
// 				ceil = points[jj+1];
// 			}
// 		}
// 	else 
// 	{
// 			color1 = numColors-3;
// 			color2 = numColors-2;
// 			color3 = numColors-1;
// 			p1 = points[numColors-3];
// 			p2 = points[numColors-2];
// 			p3 = points[numColors-1];
// 	}
	
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

	/* Quadratic color weighting */
	/*c1 = (val - p2)*(val - p3)/((p1 - p2)*(p1 - p3));
	c2 = (val - p1)*(val - p3)/((p2 - p1)*(p2 - p3));
	c3 = (val - p1)*(val - p2)/((p3 - p1)*(p3 - p2));
	
	r = c1*colors[color1][0] + c2*colors[color2][0] + c3*colors[color3][0];
	g = c1*colors[color1][1] + c2*colors[color2][1] + c3*colors[color3][1];
	b = c1*colors[color1][2] + c2*colors[color2][2] + c3*colors[color3][2];*/
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
