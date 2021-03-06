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
void GetColorForValue_MaizeBlueLatStripes(double val, double &r, double &g, double &b, double &min, double &max);
void GetColorForValue_NCLDetail(double val, double &r, double &g, double &b, double &min,double &max);
void GetColorForValue_MaizeBlue(double val, double &r, double &g, double &b, double &min, double &max);

int main( int argc, const char *argv[]){
	// Check for 3 arguments (script name, plus 2 user arguments)
	if ( argc != 3){
		cout << "Usage : "<< argv[0] << " <jobtitle | frame total>\n";
		return 1;
	};
	char vtkFileRoot[256];
	int finalFrame;
	strcpy(vtkFileRoot,argv[1]);
	istringstream ss(argv[2]);
	if ( !(ss>>finalFrame)) cerr << "ERROR: Invalid final frame number.";
	
	/* CHANGE DATA PARAMETERS */
		double dt = 0.0025;
		double relVortMin = -5.0;
		double relVortMax = 10.0;
		char scalarsDataName[16] = "relVortPanel";
		char scalarsDataTitle[64] = "Rel. Vort.";
		char scalarsDataName2[16] = "Tracer1";
		char scalarsDataTitle2[64] = "Initial Latitude";
		char titleString1[64] = "Gaussian vortices";
		char titleString2[64];
	/* END CHANGE PARAMETERS */
	
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
	
	// Setup VTK File Readers
	
	vtkPolyDataReader *sphereData[finalFrame+1];
	vtkPolyDataReader *sphereData2[finalFrame+1];
	for ( int frameJ = 0; frameJ<=finalFrame; frameJ++){
		
		j = sprintf(counterString,"%04d",frameJ);
		strcpy(vtkFileName,vtkFileRoot);
		strcat(vtkFileName,counterString);
		strcat(vtkFileName,".vtk");
		
		//initialize readers
		sphereData[frameJ]=vtkPolyDataReader::New();
		sphereData[frameJ]->SetFileName(vtkFileName);
		sphereData[frameJ]->SetScalarsName(scalarsDataName);
		sphereData2[frameJ]=vtkPolyDataReader::New();
		sphereData2[frameJ]->SetFileName(vtkFileName);
		sphereData2[frameJ]->SetScalarsName(scalarsDataName2);
	}
	cout << "... done.\n";
	
	cout << "Building color bars...\n";
	// set up vorticity colorbar
	vtkLookupTable *vortScale=vtkLookupTable::New();
	vortScale->SetScaleToLinear();
	int numColors = 1024;
	vortScale->SetRange(relVortMin,relVortMax);
	double range = relVortMax-relVortMin;
	vortScale->SetNumberOfTableValues(numColors);
	double r, g, b, val;
	for (int i=0;i<numColors;i++){
		val = relVortMin + ((double)i/numColors)*range;
		GetColorForValue_BlueWhiteRed(val,r,g, b,relVortMin,relVortMax);
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
		//GetColorForValue_JetLat(val,r,g,b,latMin,latMax);
		latScale->SetTableValue(i,r,g,b);
	}
	latScale->Build();
	
	// Setup VTK Pipeline
	vtkPolyDataMapper *sphereMapper=vtkPolyDataMapper::New();
		sphereMapper->SetScalarRange(relVortMin,relVortMax);
		sphereMapper->SetLookupTable(vortScale);
	vtkPolyDataMapper *sphereMapper2=vtkPolyDataMapper::New();
		sphereMapper2->SetScalarRange(latMin,latMax);
		sphereMapper2->SetLookupTable(latScale);
	vtkActor *sphereActor = vtkActor::New();
		sphereActor->SetMapper(sphereMapper);
	vtkActor *sphereActor2 = vtkActor::New();
		sphereActor2->SetMapper(sphereMapper2);
		sphereActor2->GetProperty()->EdgeVisibilityOn();
		sphereActor2->GetProperty()->SetEdgeColor(0,0,0);
	vtkScalarBarActor *colorbar=vtkScalarBarActor::New();
		colorbar->SetTitle(scalarsDataTitle);
		colorbar->GetTitleTextProperty()->SetColor(0.0,0.0,0.0);
		colorbar->GetLabelTextProperty()->SetColor(0.0,0.0,0.0);
		colorbar->SetOrientationToVertical();
		colorbar->SetPosition(0.025,0.1);
		colorbar->SetPosition2(0.9,0.9);
		colorbar->SetWidth(0.2);
		colorbar->SetLookupTable(vortScale);
	vtkScalarBarActor *colorbar2=vtkScalarBarActor::New();
		colorbar2->SetTitle(scalarsDataTitle2);
		colorbar2->GetTitleTextProperty()->SetColor(0.0,0.0,0.0);
		colorbar2->GetLabelTextProperty()->SetColor(0.0,0.0,0.0);
		colorbar2->SetOrientationToVertical();
		colorbar2->SetPosition(0.025,0.1);
		colorbar2->SetPosition2(0.9,0.9);
		colorbar2->SetWidth(0.2);
		colorbar2->SetLookupTable(latScale);
	vtkTextActor *textActor1 = vtkTextActor::New();
		textActor1->SetTextScaleModeToProp();
		textActor1->SetPosition(0.1,0.1);
		textActor1->SetPosition2(0.8,0.3);
		textActor1->GetTextProperty()->SetColor(0.0,0.0,0.0);
	vtkCamera *camera = vtkCamera::New();
		camera->SetPosition(5.0,0.0,1.0);
		camera->SetFocalPoint(0.0,0.0,0.0);
		camera->SetViewUp(0.0,0.0,1.0);
	vtkRenderer *renderer = vtkRenderer::New();
		renderer->AddActor(sphereActor);
		renderer->AddActor(textActor1);
		renderer->AddActor(colorbar);
		renderer->SetActiveCamera(camera);
		renderer->SetBackground(0.9,0.9,0.9);
		renderer->SetViewport(0.0,0.0,0.5,1.0);			
	vtkRenderer *renderer2 = vtkRenderer::New();
		renderer2->AddActor(sphereActor2);
		//renderer->AddActor(textActor2);
		renderer2->AddActor(colorbar2);
		renderer2->SetActiveCamera(camera);
		renderer2->SetBackground(0.9,0.9,0.9);
		renderer2->SetViewport(0.5,0.0,1.0,1.0);
	vtkRenderWindow *renWin=vtkRenderWindow::New();
		renWin->AddRenderer(renderer);
		renWin->AddRenderer(renderer2);
		renWin->SetSize(1200,600);
	vtkWindowToImageFilter *win2im = vtkWindowToImageFilter::New();
		win2im->SetInput(renWin);
	vtkJPEGWriter *writer=vtkJPEGWriter::New();
		char jpgFileRoot[128];
		char jpgFileName[128];
		strcpy(jpgFileRoot,"jpgOut/twoVorts_");
	cout << "... done.\n";
	
	
	// RENDERING LOOP
	cout << "... RENDERING ... \n";
	int nActive = 0;
	for ( int frameJ =0; frameJ<=finalFrame; frameJ++){
		
		sphereMapper->SetInput(sphereData[frameJ]->GetOutput());
		sphereMapper2->SetInput(sphereData2[frameJ]->GetOutput());
		
		j=sprintf(titleString2,"t = %7.4f, N = %6d",frameJ*dt,nActive);
		textActor1->SetInput(titleString2);
		
		j=sprintf(counterString,"%04d",frameJ);
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
	sphereMapper2->SetInput(sphereData2[finalFrame]->GetOutput());
	j=sprintf(counterString,"%04d",finalFrame+1);
	strcpy(jpgFileName,jpgFileRoot);
	strcat(jpgFileName,counterString);
	strcat(jpgFileName,".jpg");
	renWin->Render();
	win2im->Modified();
	writer->SetInput(win2im->GetOutput());
	writer->SetFileName(jpgFileName);
	writer->Write();
	
	// Clean up
	writer->Delete();
	win2im->Delete();
	renWin->Delete();
	renderer->Delete();
	renderer2->Delete();
	camera->Delete();
	sphereActor->Delete();
	sphereActor2->Delete();
	textActor1->Delete();
	colorbar->Delete();
	colorbar2->Delete();
	vortScale->Delete();
	latScale->Delete();
	sphereMapper->Delete();
	sphereMapper2->Delete();
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

void GetColorForValue_MaizeBlueLatStripes(double val, double &r, double &g, double &b, double &min,double &max)
{
	const int numColors = 3;
	const double range = max - min;
	const double colors[numColors][3] = { 0.0, 0.0, 0.4, // midnight blue
										  1.0, 1.0, 1.0, // white
										  1.0, 0.8, 0.2}; // yellow
										  
	int nStripes = 24;					  

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
		double floorj;
		double ceilj;
		double range = max - min;
		for (int j=0; j<nStripes; j++){
			floorj = min + ((double)j/(nStripes-1))*range;
			ceilj = min + ((double)(j+1)/(nStripes-1))*range;
			if (( val > floorj) && (val <= ceilj)){
				r = colors[j%numColors][0];
				g = colors[j%numColors][1];
				b = colors[j%numColors][2];
			}
		}
	}
};									  

void GetColorForValue_NCLDetail(double val, double &r, double &g, double &b, double &min,double &max)
{
	const int numColors = 254;
	const double range = max - min;
	const double colors[numColors][3] = {  0.003921568627451,   0.078431372549020,   0.996078431372549,
   0.007843137254902,   0.156862745098039,   0.992156862745098,
   0.011764705882353,   0.235294117647059,   0.988235294117647,
   0.015686274509804,   0.313725490196078,   0.984313725490196,
   0.019607843137255,   0.447058823529412,   0.980392156862745,
   0.023529411764706,   0.568627450980392,   0.976470588235294,
   0.027450980392157,   0.686274509803922,   0.972549019607843,
   0.031372549019608,   0.792156862745098,   0.968627450980392,
   0.035294117647059,   0.882352941176471,   0.964705882352941,
   0.039215686274510,   0.945098039215686,   0.960784313725490,
   0.043137254901961,   0.984313725490196,   0.956862745098039,
   0.047058823529412,   0.992156862745098,   0.952941176470588,
   0.050980392156863,   0.972549019607843,   0.949019607843137,
   0.054901960784314,   0.921568627450980,   0.945098039215686,
   0.058823529411765,   0.847058823529412,   0.941176470588235,
   0.062745098039216,   0.752941176470588,   0.937254901960784,
   0.066666666666667,   0.639215686274510,   0.933333333333333,
   0.070588235294118,   0.521568627450980,   0.929411764705882,
   0.074509803921569,   0.400000000000000,   0.925490196078431,
   0.078431372549020,   0.282352941176471,   0.921568627450980,
   0.082352941176471,   0.180392156862745,   0.917647058823529,
   0.086274509803922,   0.094117647058824,   0.913725490196078,
   0.090196078431373,   0.035294117647059,   0.909803921568627,
   0.094117647058824,   0.003921568627451,   0.905882352941176,
   0.098039215686275,                   0,   0.901960784313726,
   0.101960784313725,   0.027450980392157,   0.898039215686275,
   0.105882352941176,   0.082352941176471,   0.894117647058824,
   0.109803921568627,   0.160784313725490,   0.890196078431372,
   0.113725490196078,   0.262745098039216,   0.886274509803922,
   0.117647058823529,   0.376470588235294,   0.882352941176471,
   0.121568627450980,   0.494117647058824,   0.878431372549020,
   0.125490196078431,   0.615686274509804,   0.874509803921569,
   0.129411764705882,   0.729411764705882,   0.870588235294118,
   0.133333333333333,   0.831372549019608,   0.866666666666667,
   0.137254901960784,   0.909803921568627,   0.862745098039216,
   0.141176470588235,   0.964705882352941,   0.858823529411765,
   0.145098039215686,   0.992156862745098,   0.854901960784314,
   0.149019607843137,   0.988235294117647,   0.850980392156863,
   0.152941176470588,   0.956862745098039,   0.847058823529412,
   0.156862745098039,   0.898039215686275,   0.843137254901961,
   0.160784313725490,   0.811764705882353,   0.839215686274510,
   0.164705882352941,   0.709803921568627,   0.835294117647059,
   0.168627450980392,   0.592156862745098,   0.831372549019608,
   0.172549019607843,   0.470588235294118,   0.827450980392157,
   0.176470588235294,   0.352941176470588,   0.823529411764706,
   0.180392156862745,   0.239215686274510,   0.819607843137255,
   0.184313725490196,   0.145098039215686,   0.815686274509804,
   0.188235294117647,   0.070588235294118,   0.811764705882353,
   0.192156862745098,   0.019607843137255,   0.807843137254902,
   0.196078431372549,                   0,   0.803921568627451,
   0.200000000000000,   0.007843137254902,   0.800000000000000,
   0.203921568627451,   0.047058823529412,   0.796078431372549,
   0.207843137254902,   0.109803921568627,   0.792156862745098,
   0.211764705882353,   0.200000000000000,   0.788235294117647,
   0.215686274509804,   0.305882352941176,   0.784313725490196,
   0.219607843137255,   0.423529411764706,   0.780392156862745,
   0.223529411764706,   0.545098039215686,   0.776470588235294,
   0.227450980392157,   0.662745098039216,   0.772549019607843,
   0.231372549019608,   0.772549019607843,   0.768627450980392,
   0.235294117647059,   0.866666666666667,   0.764705882352941,
   0.239215686274510,   0.937254901960784,   0.760784313725490,
   0.243137254901961,   0.980392156862745,   0.756862745098039,
   0.247058823529412,   0.996078431372549,   0.752941176470588,
   0.250980392156863,   0.980392156862745,   0.749019607843137,
   0.254901960784314,   0.937254901960784,   0.745098039215686,
   0.258823529411765,   0.866666666666667,   0.741176470588235,
   0.262745098039216,   0.772549019607843,   0.737254901960784,
   0.266666666666667,   0.662745098039216,   0.733333333333333,
   0.270588235294118,   0.545098039215686,   0.729411764705882,
   0.274509803921569,   0.423529411764706,   0.725490196078431,
   0.278431372549020,   0.305882352941176,   0.721568627450980,
   0.282352941176471,   0.200000000000000,   0.717647058823529,
   0.286274509803922,   0.109803921568627,   0.713725490196078,
   0.290196078431373,   0.047058823529412,   0.709803921568627,
   0.294117647058824,   0.007843137254902,   0.705882352941177,
   0.298039215686275,                   0,   0.701960784313725,
   0.301960784313725,   0.019607843137255,   0.698039215686274,
   0.305882352941176,   0.070588235294118,   0.694117647058824,
   0.309803921568627,   0.145098039215686,   0.690196078431373,
   0.313725490196078,   0.239215686274510,   0.686274509803922,
   0.317647058823529,   0.352941176470588,   0.682352941176471,
   0.321568627450980,   0.470588235294118,   0.678431372549020,
   0.325490196078431,   0.592156862745098,   0.674509803921569,
   0.329411764705882,   0.709803921568627,   0.670588235294118,
   0.333333333333333,   0.811764705882353,   0.666666666666667,
   0.337254901960784,   0.898039215686275,   0.662745098039216,
   0.341176470588235,   0.956862745098039,   0.658823529411765,
   0.345098039215686,   0.988235294117647,   0.654901960784314,
   0.349019607843137,   0.992156862745098,   0.650980392156863,
   0.352941176470588,   0.964705882352941,   0.647058823529412,
   0.356862745098039,   0.909803921568627,   0.643137254901961,
   0.360784313725490,   0.831372549019608,   0.639215686274510,
   0.364705882352941,   0.729411764705882,   0.635294117647059,
   0.368627450980392,   0.615686274509804,   0.631372549019608,
   0.372549019607843,   0.498039215686275,   0.627450980392157,
   0.376470588235294,   0.376470588235294,   0.623529411764706,
   0.380392156862745,   0.262745098039216,   0.619607843137255,
   0.384313725490196,   0.160784313725490,   0.615686274509804,
   0.388235294117647,   0.082352941176471,   0.611764705882353,
   0.392156862745098,   0.027450980392157,   0.607843137254902,
   0.396078431372549,                   0,   0.603921568627451,
   0.400000000000000,   0.003921568627451,   0.600000000000000,
   0.403921568627451,   0.035294117647059,   0.596078431372549,
   0.407843137254902,   0.094117647058824,   0.592156862745098,
   0.411764705882353,   0.180392156862745,   0.588235294117647,
   0.415686274509804,   0.282352941176471,   0.584313725490196,
   0.419607843137255,   0.400000000000000,   0.580392156862745,
   0.423529411764706,   0.521568627450980,   0.576470588235294,
   0.427450980392157,   0.639215686274510,   0.572549019607843,
   0.431372549019608,   0.752941176470588,   0.568627450980392,
   0.435294117647059,   0.847058823529412,   0.564705882352941,
   0.439215686274510,   0.921568627450980,   0.560784313725490,
   0.443137254901961,   0.972549019607843,   0.556862745098039,
   0.447058823529412,   0.992156862745098,   0.552941176470588,
   0.450980392156863,   0.984313725490196,   0.549019607843137,
   0.454901960784314,   0.945098039215686,   0.545098039215686,
   0.458823529411765,   0.882352941176471,   0.541176470588235,
   0.462745098039216,   0.792156862745098,   0.537254901960784,
   0.466666666666667,   0.686274509803922,   0.533333333333333,
   0.470588235294118,   0.568627450980392,   0.529411764705882,
   0.474509803921569,   0.447058823529412,   0.525490196078431,
   0.478431372549020,   0.329411764705882,   0.521568627450980,
   0.482352941176471,   0.219607843137255,   0.517647058823529,
   0.486274509803922,   0.125490196078431,   0.513725490196078,
   0.490196078431373,   0.054901960784314,   0.509803921568627,
   0.494117647058824,   0.011764705882353,   0.505882352941176,
   0.498039215686275,                   0,   0.501960784313725,
   0.501960784313725,   0.011764705882353,   0.498039215686275,
   0.505882352941176,   0.054901960784314,   0.494117647058824,
   0.509803921568627,   0.125490196078431,   0.490196078431373,
   0.513725490196078,   0.219607843137255,   0.486274509803922,
   0.517647058823529,   0.329411764705882,   0.482352941176471,
   0.521568627450980,   0.447058823529412,   0.478431372549020,
   0.525490196078431,   0.568627450980392,   0.474509803921569,
   0.529411764705882,   0.686274509803922,   0.470588235294118,
   0.533333333333333,   0.792156862745098,   0.466666666666667,
   0.537254901960784,   0.882352941176471,   0.462745098039216,
   0.541176470588235,   0.945098039215686,   0.458823529411765,
   0.545098039215686,   0.984313725490196,   0.454901960784314,
   0.549019607843137,   0.992156862745098,   0.450980392156863,
   0.552941176470588,   0.972549019607843,   0.447058823529412,
   0.556862745098039,   0.921568627450980,   0.443137254901961,
   0.560784313725490,   0.847058823529412,   0.439215686274510,
   0.564705882352941,   0.752941176470588,   0.435294117647059,
   0.568627450980392,   0.639215686274510,   0.431372549019608,
   0.572549019607843,   0.521568627450980,   0.427450980392157,
   0.576470588235294,   0.400000000000000,   0.423529411764706,
   0.580392156862745,   0.282352941176471,   0.419607843137255,
   0.584313725490196,   0.180392156862745,   0.415686274509804,
   0.588235294117647,   0.094117647058824,   0.411764705882353,
   0.592156862745098,   0.035294117647059,   0.407843137254902,
   0.596078431372549,   0.003921568627451,   0.403921568627451,
   0.600000000000000,                   0,   0.400000000000000,
   0.603921568627451,   0.027450980392157,   0.396078431372549,
   0.607843137254902,   0.082352941176471,   0.392156862745098,
   0.611764705882353,   0.160784313725490,   0.388235294117647,
   0.615686274509804,   0.262745098039216,   0.384313725490196,
   0.619607843137255,   0.376470588235294,   0.380392156862745,
   0.623529411764706,   0.494117647058824,   0.376470588235294,
   0.627450980392157,   0.615686274509804,   0.372549019607843,
   0.631372549019608,   0.729411764705882,   0.368627450980392,
   0.635294117647059,   0.831372549019608,   0.364705882352941,
   0.639215686274510,   0.909803921568627,   0.360784313725490,
   0.643137254901961,   0.964705882352941,   0.356862745098039,
   0.647058823529412,   0.992156862745098,   0.352941176470588,
   0.650980392156863,   0.988235294117647,   0.349019607843137,
   0.654901960784314,   0.956862745098039,   0.345098039215686,
   0.658823529411765,   0.898039215686275,   0.341176470588235,
   0.662745098039216,   0.811764705882353,   0.337254901960784,
   0.666666666666667,   0.709803921568627,   0.333333333333333,
   0.670588235294118,   0.592156862745098,   0.329411764705882,
   0.674509803921569,   0.470588235294118,   0.325490196078431,
   0.678431372549020,   0.352941176470588,   0.321568627450980,
   0.682352941176471,   0.239215686274510,   0.317647058823529,
   0.686274509803922,   0.145098039215686,   0.313725490196078,
   0.690196078431373,   0.070588235294118,   0.309803921568627,
   0.694117647058824,   0.019607843137255,   0.305882352941176,
   0.698039215686274,                   0,   0.301960784313725,
   0.701960784313725,   0.007843137254902,   0.298039215686275,
   0.705882352941177,   0.047058823529412,   0.294117647058824,
   0.709803921568627,   0.109803921568627,   0.290196078431373,
   0.713725490196078,   0.200000000000000,   0.286274509803922,
   0.717647058823529,   0.305882352941176,   0.282352941176471,
   0.721568627450980,   0.423529411764706,   0.278431372549020,
   0.725490196078431,   0.545098039215686,   0.274509803921569,
   0.729411764705882,   0.662745098039216,   0.270588235294118,
   0.733333333333333,   0.772549019607843,   0.266666666666667,
   0.737254901960784,   0.866666666666667,   0.262745098039216,
   0.741176470588235,   0.937254901960784,   0.258823529411765,
   0.745098039215686,   0.980392156862745,   0.254901960784314,
   0.749019607843137,   0.996078431372549,   0.250980392156863,
   0.752941176470588,   0.980392156862745,   0.247058823529412,
   0.756862745098039,   0.937254901960784,   0.243137254901961,
   0.760784313725490,   0.866666666666667,   0.239215686274510,
   0.764705882352941,   0.772549019607843,   0.235294117647059,
   0.768627450980392,   0.662745098039216,   0.231372549019608,
   0.772549019607843,   0.545098039215686,   0.227450980392157,
   0.776470588235294,   0.423529411764706,   0.223529411764706,
   0.780392156862745,   0.305882352941176,   0.219607843137255,
   0.784313725490196,   0.200000000000000,   0.215686274509804,
   0.788235294117647,   0.109803921568627,   0.211764705882353,
   0.792156862745098,   0.047058823529412,   0.207843137254902,
   0.796078431372549,   0.007843137254902,   0.203921568627451,
   0.800000000000000,                   0,   0.200000000000000,
   0.803921568627451,   0.019607843137255,   0.196078431372549,
   0.807843137254902,   0.070588235294118,   0.192156862745098,
   0.811764705882353,   0.145098039215686,   0.188235294117647,
   0.815686274509804,   0.239215686274510,   0.184313725490196,
   0.819607843137255,   0.352941176470588,   0.180392156862745,
   0.823529411764706,   0.470588235294118,   0.176470588235294,
   0.827450980392157,   0.592156862745098,   0.172549019607843,
   0.831372549019608,   0.709803921568627,   0.168627450980392,
   0.835294117647059,   0.811764705882353,   0.164705882352941,
   0.839215686274510,   0.898039215686275,   0.160784313725490,
   0.843137254901961,   0.956862745098039,   0.156862745098039,
   0.847058823529412,   0.988235294117647,   0.152941176470588,
   0.850980392156863,   0.992156862745098,   0.149019607843137,
   0.854901960784314,   0.964705882352941,   0.145098039215686,
   0.858823529411765,   0.909803921568627,   0.141176470588235,
   0.862745098039216,   0.831372549019608,   0.137254901960784,
   0.866666666666667,   0.729411764705882,   0.133333333333333,
   0.870588235294118,   0.615686274509804,   0.129411764705882,
   0.874509803921569,   0.498039215686275,   0.125490196078431,
   0.878431372549020,   0.376470588235294,   0.121568627450980,
   0.882352941176471,   0.262745098039216,   0.117647058823529,
   0.886274509803922,   0.160784313725490,   0.113725490196078,
   0.890196078431372,   0.082352941176471,   0.109803921568627,
   0.894117647058824,   0.027450980392157,   0.105882352941176,
   0.898039215686275,                   0,   0.101960784313725,
   0.901960784313726,   0.003921568627451,   0.098039215686275,
   0.905882352941176,   0.035294117647059,   0.094117647058824,
   0.909803921568627,   0.094117647058824,   0.090196078431373,
   0.913725490196078,   0.180392156862745,   0.086274509803922,
   0.917647058823529,   0.282352941176471,   0.082352941176471,
   0.921568627450980,   0.400000000000000,   0.078431372549020,
   0.925490196078431,   0.521568627450980,   0.074509803921569,
   0.929411764705882,   0.639215686274510,   0.070588235294118,
   0.933333333333333,   0.752941176470588,   0.066666666666667,
   0.937254901960784,   0.847058823529412,   0.062745098039216,
   0.941176470588235,   0.921568627450980,   0.058823529411765,
   0.945098039215686,   0.972549019607843,   0.054901960784314,
   0.949019607843137,   0.992156862745098,   0.050980392156863,
   0.952941176470588,   0.984313725490196,   0.047058823529412,
   0.956862745098039,   0.945098039215686,   0.043137254901961,
   0.960784313725490,   0.882352941176471,   0.039215686274510,
   0.964705882352941,   0.792156862745098,   0.035294117647059,
   0.968627450980392,   0.686274509803922,   0.031372549019608,
   0.972549019607843,   0.568627450980392,   0.027450980392157,
   0.976470588235294,   0.447058823529412,   0.023529411764706,
   0.980392156862745,   0.329411764705882,   0.019607843137255,
   0.984313725490196,   0.219607843137255,   0.015686274509804,
   0.988235294117647,   0.125490196078431,   0.011764705882353,
   0.992156862745098,   0.054901960784314,   0.007843137254902,
   1.000000000000000,                   0,   0.003921568627451};

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
		double floorj;
		double ceilj;
		double range = max - min;
		for (int j=0; j<numColors; j++){
			floorj = min + ((double)j/(numColors-1))*range;
			ceilj = min + ((double)(j+1)/(numColors-1))*range;
			if (( val > floorj) && (val <= ceilj)){
				r = colors[j][0];
				g = colors[j][1];
				b = colors[j][2];
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
	//double points[9] = {-1.5708, -0.7854, 0.0, 0.6, 0.7, 0.7854, 0.9, 1.0, 1.5708};
	//					cyan      blue    orchid l.sl.blue med.turq green  gold	 dk. org. yellow
	double points[9] = {-1.5708, -0.7854, 0.0,   0.17,     0.265,   0.35,  0.9,  1.0,     1.5708};
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
