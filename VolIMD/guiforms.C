/*                         
                           (C) 1997
              Computer Centre University of Stuttgart
                         Allmandring 30
                       D-70550 Stuttgart
                            Germany

THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

Permission to use, copy, modify and distribute this software is hereby 
granted without fee, provided that the above and this permission notice 
appear in all copies of this software and that the authorship is 
acknowledged.

*/
/* program name ..............  VolRend  
   file name .................  guiforms.C
   author ....................  Roland Niemeier
   version ...................  0.9beta0
   date of last change .......  07/25/97
   description : interface for direct volume visualization of scalar fields. 
                 The application development is intended for 
                 the interactive visualization of volume data sets. 
                 It is based on OpenInventor, VolPack, Xforms.
   purpose of guiforms.C: Collection of xforms functions for popup menus
*/

#include <forms.h>
#include <volpack.h>
#include "global.h"
#include "render.h"
#include <stdlib.h>
#include <math.h>

// declaration of global variables
/* declaration of number of ramp points for building classification ramps */
extern int n_scalar, n_gradient;
/* arrays for storage for ramp point values for building classification ramps */
extern float SRX[256], SRY[256], GRX[256], GRY[256];
/* it is just for VP_SCALAR_MAX that we include volpack.h */
/* arrays for storage for ramp point values for building weight ramps */
extern float WX[VP_SCALAR_MAX+1], WY[VP_SCALAR_MAX+1];
/* name of server where the simulation is running */
extern char server_name[256];
/* for hints in the histogram display the maxima and minima of a received
   data set is stored in last_received... */
extern float last_received_kinetic_minimum, last_received_kinetic_maximum;
extern float last_received_potential_minimum, last_received_potential_maximum;
/* boolean for determining if kinetic or potential energy data 
   is displayed used here for histogramm display */
extern unsigned char KineticEnergyFile;
/* Boolean variable to force redraws */
extern unsigned char FORCEDRAW;
extern unsigned short base_port;

// definition of local functions
void show_message(const char *s1, const char *s2, const char *s3);
// definition of file local variables
// filename for storing or loading material properties 
char material_filename[256];

// There are two different ways to process events with xforms:
//  1. events are processed all in an internal loop of a subroutine.
//  2. events are processed by callbacks. For this structures must be defined
//     outside the subroutine.
//  In general I prefere the second possibility (not allocating memory for
//  variables that are not used) but one example for the first way is the
//  classification functions menu for new classifications.
//  

// Structure used by Callback Routines called from active_xyplot for the 
// classification functions (assignmment of opacities)  
typedef struct {
        FL_FORM *Classification;
        FL_OBJECT *ScalarClassificationPlot;
        FL_OBJECT *GradientClassificationPlot;
        FL_OBJECT *DisplayValues;
        FL_OBJECT *ApplyBox;
        FL_OBJECT *SaveBox;
        FL_OBJECT *LoadBox;
        FL_OBJECT *QuitBox;
        FL_OBJECT *status;
	FL_OBJECT *GradientPointsCounter;
	FL_OBJECT *ScalarPointsCounter;
	FL_OBJECT *ClassificationTextBox;
	FL_OBJECT *InformationTextBox1;
	FL_OBJECT *InformationTextBox2;
        void *vdata;
        long ldata;
} FD_Classification;
FD_Classification *create_form_Classification(void);
FD_Classification *xypui;

typedef struct{
  FL_FORM *histogram_form;
  FL_OBJECT *histogram_xyplot;
  FL_OBJECT *MinimumHint;
  FL_OBJECT *MaximumHint;
  //FL_OBJECT *histogram_quit;  
} FD_histogram;
FD_histogram *histogram_object;

void forms_initialization(int *argc, char **argv){
/************************************************************************/
/* Makes some initial settings for using xforms				*/
/************************************************************************/
   printf("forms_initialization ...\n");
   fl_initialize(argc, argv, "dummy", NULL, 0);
}

/*****************************************************************************/
/* Callback Routines called from routines defined in this file		     */
/*****************************************************************************/

void value_cb(FL_OBJECT *ob, long)
/*****************************************************************************/
/* Print the x-y-values of the selected point (called from active_xyplot)    */
/*****************************************************************************/
{
    float x, y;
    int i;
    char buf[64];
    fl_get_xyplot(ob, &x, &y, &i);
    if (i < 0) return;
    sprintf(buf,"X=%.2f  Y=%.2f",x,y);
    fl_set_object_label(xypui->DisplayValues, buf);
}

void SetScalarPointsNumberCB(FL_OBJECT *ob, long)
/*****************************************************************************/
/* for changing the number of points in the plot that define the     	     */
/* function scalar_value --> opacity_1 (called from active_xyplot)           */
/*****************************************************************************/
{
int number_of_points, new_number_of_points, i, counter=0;
float maximum_distance, minimum_distance, compare_distance;
float SRX_counter, SRY_counter;
// get the new number of points
  new_number_of_points = (int) (fl_get_counter_value(ob)+0.5);
  printf("SetScalarPointsNumberCB: Number of Points is %d\n",new_number_of_points);
  fl_get_xyplot_data(xypui->ScalarClassificationPlot, SRX, SRY, &number_of_points);
  if (new_number_of_points > number_of_points){
    //insert points between points that have maximum difference in their x-coordinates
    while (new_number_of_points > number_of_points){
      // determine maximum distance
      maximum_distance = SRX[1] - SRX[0];
      for(i=1;i<number_of_points;i++){
        compare_distance = SRX[i+1] - SRX[i];
        if (maximum_distance < compare_distance){
          counter=i;
          maximum_distance = compare_distance;
        }
      }
      // insert points
      SRX_counter = 0.5 * (SRX[counter+1] + SRX[counter]); 
      SRY_counter = 0.5 * (SRY[counter+1] + SRY[counter]); 
      number_of_points++;
      for (i=number_of_points; i>(counter+1); i--){
        SRX[i] = SRX[i-1]; SRY[i] = SRY[i-1];
      }
      SRX[counter+1] = SRX_counter; SRY[counter+1] = SRY_counter;
    }    
  }
  else if ((new_number_of_points < number_of_points) && (number_of_points > 2)){
    //determine minimum distance between the neighbouring points of a point
    while (new_number_of_points < number_of_points){
      minimum_distance = SRX[2] - SRX[0]; counter = 1;
      for(i=2;i<(number_of_points-1);i++){
        compare_distance = SRX[i+1] - SRX[i-1];
        if (minimum_distance > compare_distance){
          counter=i;
          minimum_distance = compare_distance;
        }
      }
    //delete points
      number_of_points--;
      for (i=counter; i<number_of_points; i++){
        SRX[i] = SRX[i+1]; SRY[i] = SRY[i+1];
      }
    }    
  }
  for (i=0;i<number_of_points;i++) {
    printf("Scalar Ramp Point %d X Value: %.2f\n",i+1,SRX[i]);
    printf("Scalar Ramp Point %d Y Value: %.2f\n",i+1,SRY[i]);
  }
  fl_set_xyplot_data(xypui->ScalarClassificationPlot, SRX, SRY, number_of_points, "","","");
  n_scalar = number_of_points;
}

void SetGradientPointsNumberCB(FL_OBJECT *ob, long)
/*****************************************************************************/
/* for changing the number of points in the plot that define the     	     */
/* function gradient_value --> opacity_2 (called from active_xyplot)         */
/*****************************************************************************/
{
int number_of_points, new_number_of_points, i, counter=0;
float maximum_distance, minimum_distance, compare_distance;
float GRX_counter, GRY_counter;
// get the new number of points
  new_number_of_points = (int) (fl_get_counter_value(ob)+0.5);
  fl_get_xyplot_data(xypui->GradientClassificationPlot, GRX, GRY, &number_of_points);
  if (new_number_of_points > number_of_points){
    //insert points between points that have maximum distance
    while (new_number_of_points > number_of_points){
      // determine maximum distance
      maximum_distance = GRX[1] - GRX[0];
      for(i=1;i<number_of_points;i++){
        compare_distance = GRX[i+1] - GRX[i];
        if (maximum_distance < compare_distance){
          counter=i;
          maximum_distance = compare_distance;
        }
      }
      // insert points
      GRX_counter = 0.5 * (GRX[counter+1] + GRX[counter]); 
      GRY_counter = 0.5 * (GRY[counter+1] + GRY[counter]); 
      number_of_points++;
      for (i=number_of_points; i>(counter+1); i--){
        GRX[i] = GRX[i-1]; GRY[i] = GRY[i-1];
      }
      GRX[counter+1] = GRX_counter; GRY[counter+1] = GRY_counter;
    }    
  }
  else if ((new_number_of_points < number_of_points) && (number_of_points > 2)){
    //determine minimum distance between the neighbouring points of a point
    while (new_number_of_points < number_of_points){
      minimum_distance = GRX[2] - GRX[0]; counter = 1;
      for(i=2;i<(number_of_points-1);i++){
        compare_distance = GRX[i+1] - GRX[i-1];
        if (minimum_distance > compare_distance){
          counter=i;
          minimum_distance = compare_distance;
        }
      }
    //delete points
      number_of_points--;
      for (i=counter; i<number_of_points; i++){
        GRX[i] = GRX[i+1]; GRY[i] = GRY[i+1];
      }
    }    
  }
  for (i=0;i<number_of_points;i++) {
    printf("Gradient Ramp Point %d X Value: %.2f\n",i+1,GRX[i]);
    printf("Gradient Ramp Point %d Y Value: %.2f\n",i+1,GRY[i]);
  }
  fl_set_xyplot_data(xypui->GradientClassificationPlot, GRX, GRY, number_of_points, "","","");
  n_gradient = number_of_points;
}

void SetWeightPointsNumberCB(FL_OBJECT *, int material_weight_points, int new_material_weight_points)
/*****************************************************************************/
/* for changing the number of points in the plot that define the     	     */
/* function scalar_value --> material_proportions of current_material	     */
/* (called from modify_material_popup)				             */
/*****************************************************************************/
{
int i, counter=0;
float maximum_distance, minimum_distance, compare_distance;
float WX_counter, WY_counter;
// get the new number of points
  //new_material_weight_points = (int) (fl_get_counter_value(ob)+0.5);
  if (new_material_weight_points > material_weight_points){
    //insert points between points that have maximum distance
    while (new_material_weight_points > material_weight_points){
      // determine maximum distance
      maximum_distance = WX[1] - WX[0];
      for(i=1;i<material_weight_points;i++){
        compare_distance = WX[i+1] - WX[i];
        if (maximum_distance < compare_distance){
          counter=i;
          maximum_distance = compare_distance;
        }
      }
      // insert points
      WX_counter = 0.5 * (WX[counter+1] + WX[counter]); 
      WY_counter = 0.5 * (WY[counter+1] + WY[counter]); 
      material_weight_points++;
      for (i=material_weight_points; i>(counter+1); i--){
        WX[i] = WX[i-1]; WY[i] = WY[i-1];
      }
      WX[counter+1] = WX_counter; WY[counter+1] = WY_counter;
    }    
  }
  else if ((new_material_weight_points < material_weight_points) && (material_weight_points > 2)){
    //determine minimum distance between the neighbouring points of a point
    while (new_material_weight_points < material_weight_points){
      minimum_distance = WX[2] - WX[0]; counter = 1;
      for(i=2;i<(material_weight_points-1);i++){
        compare_distance = WX[i+1] - WX[i-1];
        if (minimum_distance > compare_distance){
          counter=i;
          minimum_distance = compare_distance;
        }
      }
    //delete points
      material_weight_points--;
      for (i=counter; i<material_weight_points; i++){
        WX[i] = WX[i+1]; WY[i] = WY[i+1];
      }
    }    
  }
  //for (i=0;i<material_weight_points;i++) {
  //  printf("Weight Ramp Point %d X Value: %.2f\n",i+1,WX[i]);
  //  printf("Weight Ramp Point %d Y Value: %.2f\n",i+1,WY[i]);
  //}
}

void save_classification_file(const char *filename){
  FILE *classification_fp;
  int i, number_of_scalar_points, number_of_gradient_points;
  if ((classification_fp = fopen(filename, "w")) == NULL) {
    show_message("cannot open classification ramp file ","","");
  }
  else{
    // write scalar ramp points
    fl_get_xyplot_data(xypui->ScalarClassificationPlot, SRX, SRY, &number_of_scalar_points);
    fprintf(classification_fp,"%d\n",number_of_scalar_points);
    for(i=0;i<number_of_scalar_points;i++){
      fprintf(classification_fp,"%.2f %.2f\n",SRX[i],SRY[i]);
    }
    // write gradient ramp points
    fl_get_xyplot_data(xypui->GradientClassificationPlot, GRX, GRY, &number_of_gradient_points);
    fprintf(classification_fp,"%d\n",number_of_gradient_points);
    for(i=0;i<number_of_gradient_points;i++){
      fprintf(classification_fp,"%.2f %.2f\n",GRX[i],GRY[i]);
    }
    fclose(classification_fp);
  }
}

void load_classification_file(const char *filename){
  FILE *classification_fp;
  int i, number_of_points;
  if ((classification_fp = fopen(filename, "r")) == NULL) {
    show_message("cannot open classification ramp file ", "","");
  }
  else{
    // read scalar ramp points
    fscanf(classification_fp,"%d",&number_of_points);
    //printf("1. number of points is %d\n",number_of_points);
    for(i=0;i<number_of_points;i++){
      fscanf(classification_fp,"%f %f",&(SRX[i]),&(SRY[i]));
      //printf("%.2f %.2f\n",SRX[i],SRY[i]);
    }
    fl_set_xyplot_data(xypui->ScalarClassificationPlot, SRX, SRY, number_of_points, "","","");
    n_scalar = number_of_points;
    fl_set_counter_value(xypui->ScalarPointsCounter,(double)n_scalar);
    // read gradient ramp points
    fscanf(classification_fp,"%d",&number_of_points);
    //printf("2. number of points is %d\n",number_of_points);
    for(i=0;i<number_of_points;i++){
      fscanf(classification_fp,"%f %f",&(GRX[i]),&(GRY[i]));
      //printf("%.2f %.2f\n",GRX[i],GRY[i]);
    }
    fl_set_xyplot_data(xypui->GradientClassificationPlot, GRX, GRY, number_of_points, "","","");
    n_gradient = number_of_points;
    fl_set_counter_value(xypui->GradientPointsCounter,(double)n_gradient);
    fclose(classification_fp);
  }
}

void save_classification_ramps_cb(FL_OBJECT *, long){
  const char *filename;
  fl_set_fselector_placement(FL_PLACE_FREE);
  filename = fl_show_fselector("Save Classification File", 0, "*.cla",0);
  if (filename != NULL){
    printf("selected classification filename is %s\n",filename);
    save_classification_file(filename);
  }
  else{
    show_message("No file name selected \n","","");
  }
}

void load_classification_ramps_cb(FL_OBJECT *, long){
  const char *filename;
  fl_set_fselector_placement(FL_PLACE_FREE);
  filename = fl_show_fselector("Load Classification File", 0, "*.cla",0);
  if (filename != NULL){
    printf("selected classification filename is %s\n",filename);
    load_classification_file(filename);
  }
  else{
    show_message("No file name selected \n","","");
  }
}


void apply_cb(FL_OBJECT *, long)
/************************************************************************/
/* Get the new Ramp and Start new classification module (in graph.C)    */
/* (called from active_xyplot)						*/
/************************************************************************/
{
    int i; 
    fl_get_xyplot_data(xypui->ScalarClassificationPlot, 
           SRX, SRY, &n_scalar);
    if(n_scalar <= 0){
      printf("Can not read scalar xyplot data points\n");
      return ;
    }
    for (i=0;i<n_scalar;i++){
      printf("Scalar Point %d: X=%f  Y=%f\n",i,SRX[i],SRY[i]);
    }
    fl_get_xyplot_data(xypui->GradientClassificationPlot, 
           GRX, GRY, &n_gradient);
    if(n_gradient <= 0){
      printf("Can not read gradient xyplot data points\n");
      return ;
    }
    for (i=0;i<n_gradient;i++){
      printf("Gradient Point %d: X=%f  Y=%f\n",i,GRX[i],GRY[i]);
    }
    new_classification();
}

void active_xyplot(void)
/************************************************************************/
/* Called after pulldown selection classification from menu modify	*/
/* Manages the form							*/
/************************************************************************/
{
   xypui = create_form_Classification();

   /* fill-in form initialization code */
   fl_set_object_dblbuffer(xypui->DisplayValues, 1);

   fl_set_xyplot_data(xypui->ScalarClassificationPlot, SRX, SRY, n_scalar, "","","");
   fl_set_xyplot_xtics(xypui->ScalarClassificationPlot, 5, 2);
   fl_set_xyplot_ytics(xypui->ScalarClassificationPlot, 5, 2);
   fl_set_xyplot_ybounds(xypui->ScalarClassificationPlot, 0.0, 1.0);

   fl_set_xyplot_data(xypui->GradientClassificationPlot, GRX, GRY, n_gradient, "","","");
   fl_set_xyplot_xtics(xypui->GradientClassificationPlot, 5, 2);
   fl_set_xyplot_ytics(xypui->GradientClassificationPlot, 5, 2);
   fl_set_xyplot_ybounds(xypui->GradientClassificationPlot, 0.0, 1.0);

   /* show the first form */
    fl_show_form(xypui->Classification, FL_PLACE_MOUSE|FL_FREE_SIZE, FL_FULLBORDER, "Classification");
   /*fl_show_form(xypui->Classification,FL_PLACE_MOUSE,FL_TRANSIENT,"Classification");*/
   fl_do_forms();
   fl_hide_form(xypui->Classification);
   fl_free_form(xypui->Classification);
}


FD_Classification *create_form_Classification(void)
/************************************************************************/
/* Creates the classification form and determines the layout	        */
/************************************************************************/
{
  FL_OBJECT *obj;
  FD_Classification *fdui = (FD_Classification *) fl_calloc(1, sizeof(FD_Classification));

  fdui->Classification = fl_bgn_form(FL_NO_BOX, 480, 650);
  obj = fl_add_box(FL_UP_BOX,0,0,480,650,"");
  fdui->ScalarClassificationPlot = obj = fl_add_xyplot(FL_ACTIVE_XYPLOT,20,120,330,211,"Scalar Value Classification");
    fl_set_object_boxtype(obj,FL_FRAME_BOX);
    fl_set_object_callback(obj,value_cb,0);
  fdui->GradientClassificationPlot = obj = fl_add_xyplot(FL_ACTIVE_XYPLOT,20,361,330,211,"Scalar Gradient Classification");
    fl_set_object_boxtype(obj,FL_FRAME_BOX);
    fl_set_object_callback(obj,value_cb,0);
  fdui->DisplayValues = obj = fl_add_box(FL_BORDER_BOX,342,577,108,24,"");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->ApplyBox = obj = fl_add_button(FL_NORMAL_BUTTON,36,611,83,29,"Apply");
    fl_set_object_callback(obj,apply_cb,0);
  fdui->SaveBox = obj = fl_add_button(FL_NORMAL_BUTTON,135,611,83,29,"Save");
    fl_set_object_callback(obj,save_classification_ramps_cb,0);
  fdui->LoadBox = obj = fl_add_button(FL_NORMAL_BUTTON,235,611,83,29,"Load");
    fl_set_object_callback(obj,load_classification_ramps_cb,0);
  fdui->QuitBox = obj = fl_add_button(FL_NORMAL_BUTTON,335,611,83,29,"Quit");
  fdui->GradientPointsCounter = obj = fl_add_counter(FL_SIMPLE_COUNTER,359,505,91,29,"Number of points");
    fl_set_object_callback(obj,SetGradientPointsNumberCB,0);
    fl_set_counter_precision(obj,0);
    fl_set_counter_bounds(obj,2.0,255.0);
    fl_set_counter_value(obj,(double)n_gradient);
  fdui->ScalarPointsCounter = obj = fl_add_counter(FL_SIMPLE_COUNTER,359,264,91,29,"Number of points");
    fl_set_object_callback(obj,SetScalarPointsNumberCB,0);
    fl_set_counter_precision(obj,0);
    fl_set_counter_bounds(obj,2.0,255.0);
    fl_set_counter_value(obj,(double)n_scalar);
  fdui->ClassificationTextBox = obj = fl_add_text(FL_NORMAL_TEXT,150,10,200,40,"Classification");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  fdui->InformationTextBox1 = obj = fl_add_text(FL_NORMAL_TEXT,20,60,440,20,"sets the opacity of a  voxel  dependent on the  voxels scalar value and its gradient");
  fdui->InformationTextBox2 = obj = fl_add_text(FL_NORMAL_TEXT,20,90,200,20,"by multiplying the two y-axis values");
  fl_end_form();

  return fdui;
}


void rotation_slider(int *number_of_rotations, double *rotation_increment)
/********************************************************************************/
/* Activates the rotation values pulldown selection from the menu modify    	*/
/* for changing the angle increment of a rotation and the number of rotations	*/		
/* carried out if pulldown menu rotate around ... is used			*/
/********************************************************************************/
{
  typedef struct {
	FL_FORM *rotation_values_form;
	FL_OBJECT *RotationAngleSlider;
        FL_OBJECT *RotationNumberSlider;
        FL_OBJECT *RotationApply;
        FL_OBJECT *RotationReset;
        FL_OBJECT *RotationQuit;
	void *vdata;
	long ldata;
  } FD_rotation_values;
  double old_rotation_increment = *rotation_increment;
  int old_number_of_rotations = *number_of_rotations;
   //rotationui = create_form_rotation_slider();
  FL_OBJECT *obj;
  FD_rotation_values *fdui = (FD_rotation_values *) fl_calloc(1, sizeof(FD_rotation_values));
  fdui->rotation_values_form = fl_bgn_form(FL_NO_BOX, 300, 230);
  obj = fl_add_box(FL_UP_BOX,0,0,300,230,"");
  fdui->RotationAngleSlider = obj = fl_add_valslider(FL_HOR_SLIDER,30,30,230,30,"Rotation angle increment");
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.0,90.0);
    fl_set_slider_step(obj,0.5);
    fl_set_slider_precision(obj,1);
   fl_set_slider_value(obj,*rotation_increment);
   fdui->RotationNumberSlider = obj = fl_add_valslider(FL_HOR_SLIDER,30,90,230,30,"Number of rotations");
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.0,50.0);
    fl_set_slider_step(obj,1.0);
    fl_set_slider_precision(obj,0);
    fl_set_slider_value(obj,(double) (*number_of_rotations));
  fdui->RotationApply = obj = fl_add_button(FL_NORMAL_BUTTON,30,160,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->RotationReset = obj = fl_add_button(FL_NORMAL_BUTTON,110,160,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->RotationQuit = obj = fl_add_button(FL_NORMAL_BUTTON,190,160,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fl_end_form();
  obj = NULL;
  /* show the form */
  fl_show_form(fdui->rotation_values_form,FL_PLACE_MOUSE,FL_TRANSIENT,"Rotation");
  while (obj != fdui->RotationQuit){
    obj = fl_do_forms();
    if (obj == fdui->RotationApply){
      *number_of_rotations = (int) fl_get_slider_value(fdui->RotationNumberSlider);
      *rotation_increment =  (float)fl_get_slider_value(fdui->RotationAngleSlider);
    } 
    if (obj == fdui->RotationReset){
      fl_set_slider_value(fdui->RotationNumberSlider,(double)old_number_of_rotations);
      *number_of_rotations = old_number_of_rotations; 
      fl_set_slider_value(fdui->RotationAngleSlider,old_rotation_increment);
      *rotation_increment = old_rotation_increment;
    }
  }
  fl_hide_form(fdui->rotation_values_form);
  fl_free_form(fdui->rotation_values_form);
}


void histogram_slider(int *histogram_minimum, int *histogram_maximum)
/********************************************************************************/
/* Activates the histogram values pulldown selection from the menu modify    	*/
/* for changing the histogram boundaries (y-axis)				*/		
/* carried out if pulldown menu histogram ... is used				*/
/********************************************************************************/
{
  typedef struct {
	FL_FORM *histogram_values_form;
	FL_OBJECT *MinimumSlider;
        FL_OBJECT *MaximumSlider;
        FL_OBJECT *HistogramApply;
        FL_OBJECT *HistogramReset;
        FL_OBJECT *HistogramQuit;
	//void *vdata;
	//long ldata;
  } FD_histogram_values;
  double old_histogram_minimum = *histogram_minimum;
  int old_histogram_maximum = *histogram_maximum;
  FL_OBJECT *obj;
  FD_histogram_values *fdui = (FD_histogram_values *) fl_calloc(1, sizeof(FD_histogram_values));
  fdui->histogram_values_form = fl_bgn_form(FL_NO_BOX, 300, 230);
  obj = fl_add_box(FL_UP_BOX,0,0,300,230,"");
  fdui->MaximumSlider = obj = fl_add_valslider(FL_HOR_SLIDER,30,30,230,30,"upper boundary");
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.0,1000.0);
    fl_set_slider_step(obj,1.0);
    fl_set_slider_precision(obj,0);
    fl_set_slider_value(obj,(double) (*histogram_maximum));
  fdui->MinimumSlider = obj = fl_add_valslider(FL_HOR_SLIDER,30,90,230,30,"lower boundary");
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.0,50.0);
    fl_set_slider_step(obj,1.0);
    fl_set_slider_precision(obj,0);
    fl_set_slider_value(obj,(double) (*histogram_minimum));
  fdui->HistogramApply = obj = fl_add_button(FL_NORMAL_BUTTON,30,160,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->HistogramReset = obj = fl_add_button(FL_NORMAL_BUTTON,110,160,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->HistogramQuit = obj = fl_add_button(FL_NORMAL_BUTTON,190,160,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fl_end_form();
  obj = NULL;
  /* show the form */
  fl_show_form(fdui->histogram_values_form,FL_PLACE_MOUSE,FL_TRANSIENT,"Histogram Boundaries");
  while (obj != fdui->HistogramQuit){
    obj = fl_do_forms();
    if (obj == fdui->HistogramApply){
      *histogram_minimum = rint(fl_get_slider_value(fdui->MinimumSlider));
      *histogram_maximum = rint(fl_get_slider_value(fdui->MaximumSlider));
    } 
    if (obj == fdui->HistogramReset){
      fl_set_slider_value(fdui->MinimumSlider,(double)old_histogram_minimum);
      *histogram_minimum = old_histogram_minimum; 
      fl_set_slider_value(fdui->MaximumSlider,(double)old_histogram_maximum);
      *histogram_maximum = old_histogram_maximum;
    }
  }
  fl_hide_form(fdui->histogram_values_form);
  fl_free_form(fdui->histogram_values_form);
  if ((*histogram_minimum) > (*histogram_maximum)){
    printf("minimum and maximum %d %d\n",*histogram_minimum,*histogram_maximum);
    show_message("sorry!","minimum greater than maximum","reset to previous values");
    *histogram_minimum = old_histogram_minimum; 
    *histogram_maximum = old_histogram_maximum;
  }
}


void modify_lighting_popup(double *r_light, double *g_light, double *b_light,
                           double *x_direction, double *y_direction, 
                           double *z_direction, int current_light)
/********************************************************************************/
/* Activates the light properties pulldown selection from the menu modify    	*/
/* for changing the direction and colors of the directional light sources	*/		
/********************************************************************************/
{
typedef struct {
        FL_FORM *modify_light;
        FL_OBJECT *LightPropertiesTextBox;
        FL_OBJECT *LightTextBox;
        FL_OBJECT *ColorTextBox;
        FL_OBJECT *DirectionTextBox;
        FL_OBJECT *ColorRedSlider;
        FL_OBJECT *ColorGreenSlider;
        FL_OBJECT *ColorBlueSlider;
        FL_OBJECT *XDirectionSlider;
        FL_OBJECT *YDirectionSlider;
        FL_OBJECT *ZDirectionSlider;
        FL_OBJECT *ModifyLightApply;
        FL_OBJECT *ModifyLightReset;
        FL_OBJECT *ModifyLightQuit;
        void *vdata;
        long ldata;
} FD_light_sliders;
  FL_OBJECT *obj;
  FD_light_sliders *fdui = (FD_light_sliders *) fl_calloc(1, sizeof(FD_light_sliders));
  double 
    oldr_light = *r_light, oldg_light = *g_light, oldb_light = *b_light,
    oldx_direction = *x_direction, oldy_direction = *y_direction, 
    oldz_direction = *z_direction;
  char which_light[20];
  fdui->modify_light = fl_bgn_form(FL_NO_BOX, 650, 310);
  obj = fl_add_box(FL_UP_BOX,0,0,650, 310,"");
  fdui->LightPropertiesTextBox = obj = fl_add_text(FL_NORMAL_TEXT,220,5,200,40,"Light Properties");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  sprintf(which_light,"of Light %d",current_light);
  fdui->LightTextBox = obj = fl_add_text(FL_NORMAL_TEXT,260,38,200,40,which_light);
    fl_set_object_lsize(obj,FL_MEDIUM_SIZE);
  fdui->ColorTextBox = obj = fl_add_text(FL_NORMAL_TEXT,40,75,200,30,"Colors of light source");
    fl_set_object_boxtype(obj,FL_SHADOW_BOX);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->DirectionTextBox = obj = fl_add_text(FL_NORMAL_TEXT,320,75,170,30,"Direction of light source");
    fl_set_object_boxtype(obj,FL_SHADOW_BOX);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->ColorRedSlider = obj = fl_add_valslider(FL_HOR_SLIDER,30,125,230,30,"Red");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_RED,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetAmbientRed,0);
  fdui->ColorGreenSlider = obj = fl_add_valslider(FL_HOR_SLIDER,30,185,230,30,"Green");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_GREEN,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetAmbientGreen,0);
  fdui->ColorBlueSlider = obj = fl_add_valslider(FL_HOR_SLIDER,30,245,230,30,"Blue");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_BLUE,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetAmbientBlue,0);
  fdui->XDirectionSlider = obj = fl_add_valslider(FL_HOR_SLIDER,290,125,230,30,"x");
    fl_set_slider_bounds(obj,-1.0,1.0);
    fl_set_slider_return(obj, 0);
   // fl_set_object_color(obj,FL_RED,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetDiffuseRed,0);
  fdui->YDirectionSlider = obj = fl_add_valslider(FL_HOR_SLIDER,290,185,230,30,"y");
    fl_set_slider_bounds(obj,-1.0,1.0);
    //fl_set_slider_step(obj,1.0);
    fl_set_slider_return(obj, 0);
   // fl_set_object_color(obj,FL_GREEN,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetDiffuseGreen,0);
  fdui->ZDirectionSlider = obj = fl_add_valslider(FL_HOR_SLIDER,290,245,230,30,"z");
    fl_set_slider_bounds(obj,-1.0,1.0);
    fl_set_slider_return(obj, 0);
   // fl_set_object_color(obj,FL_BLUE,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetDiffuseBlue,0);
  fdui->ModifyLightReset = obj = fl_add_button(FL_NORMAL_BUTTON,550,185,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->ModifyLightApply = obj = fl_add_button(FL_NORMAL_BUTTON,550,125,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->ModifyLightQuit = obj = fl_add_button(FL_NORMAL_BUTTON,550,245,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  obj = fdui->ZDirectionSlider;
  fl_end_form();
  fl_set_slider_value(fdui->ColorRedSlider,oldr_light);
  fl_set_slider_value(fdui->ColorGreenSlider,oldg_light);
  fl_set_slider_value(fdui->ColorBlueSlider,oldb_light);
  fl_set_slider_value(fdui->XDirectionSlider,oldx_direction);
  fl_set_slider_value(fdui->YDirectionSlider,oldy_direction);
  fl_set_slider_value(fdui->ZDirectionSlider,oldz_direction);
  /* show the form */
  fl_show_form(fdui->modify_light,FL_PLACE_MOUSE,FL_TRANSIENT,"Modify Light");
  while (obj != fdui->ModifyLightQuit){
    obj = fl_do_forms();
    if (obj == fdui->ModifyLightApply){
      *r_light = fl_get_slider_value(fdui->ColorRedSlider); 
      *g_light = fl_get_slider_value(fdui->ColorGreenSlider);
      *b_light = fl_get_slider_value(fdui->ColorBlueSlider);
      *x_direction = fl_get_slider_value(fdui->XDirectionSlider); 
      *y_direction = fl_get_slider_value(fdui->YDirectionSlider);
      *z_direction = fl_get_slider_value(fdui->ZDirectionSlider);
      setLightColors(*r_light,*g_light,*b_light);
      setLightDirection(*x_direction,*y_direction,*z_direction);
      multi_rotate_image(1,0.0,1);
    }
    if (obj == fdui->ModifyLightReset){
      fl_set_slider_value(fdui->ColorRedSlider,oldr_light); *r_light = oldr_light;
      fl_set_slider_value(fdui->ColorGreenSlider,oldg_light); *g_light = oldg_light;
      fl_set_slider_value(fdui->ColorBlueSlider,oldb_light); *b_light = oldb_light;
      fl_set_slider_value(fdui->XDirectionSlider,oldx_direction); *x_direction = oldx_direction;
      fl_set_slider_value(fdui->YDirectionSlider,oldy_direction); *y_direction = oldy_direction;
      fl_set_slider_value(fdui->ZDirectionSlider,oldz_direction); *z_direction = oldz_direction;
      setLightColors(*r_light,*g_light,*b_light);
      setLightDirection(*x_direction,*y_direction,*z_direction);
      multi_rotate_image(1,0.0,1);
    }
  }
  fl_hide_form(fdui->modify_light);
  fl_free_form(fdui->modify_light);
}
//int save_material_properties(const char *filename){
//  printf("Selected File Name is %s\n",filename);
//  sprintf(material_filename,"%s",filename); return 1;
//}
//int load_material_properties(const char *filename){
//  printf("Selected File Name is %s\n",filename);
//  sprintf(material_filename,"%s",filename); return 1;
//}
void modify_material_popup(double *r_ambient, double *g_ambient, double *b_ambient,
		       double *r_diffuse,double *g_diffuse, double *b_diffuse,
                       double *r_specular,double *g_specular, double *b_specular,
		       double *shinyness, 
                       int *current_material, int *number_of_materials,
                       int material_weight_points)
/********************************************************************************/
/* Activates the material properties pulldown selection from the menu modify    */
/* for changing the materials ambient, diffuse and specular coefficients	*/		
/* and for interactive modification of the function: scalar value of the voxel	*/
/* -> probability that the selected material is contained in this voxel		*/
/********************************************************************************/
{
const char *filename;
typedef struct {
        FL_FORM *modify_material_form;
        FL_OBJECT *AmbientRedSlider;
        FL_OBJECT *ShinynessSlider;
        FL_OBJECT *ShinynessTextBox;
        FL_OBJECT *AmbientTextBox;
        FL_OBJECT *DiffuseTextBox;
        FL_OBJECT *SpecularTextBox;
        FL_OBJECT *MaterialPropertiesTextBox;
        FL_OBJECT *AmbientGreenSlider;
        FL_OBJECT *AmbientBlueSlider;
        FL_OBJECT *DiffuseRedSlider;
        FL_OBJECT *DiffuseGreenSlider;
        FL_OBJECT *DiffuseBlueSlider;
        FL_OBJECT *SpecularRedSlider;
        FL_OBJECT *SpecularGreenSlider;
        FL_OBJECT *SpecularBlueSlider;
        FL_OBJECT *ModifyMaterialApply;
        //FL_OBJECT *ModifyMaterialReset;
        FL_OBJECT *ModifyMaterialQuit;
        FL_OBJECT *ModifyMaterialLoad;
        FL_OBJECT *ModifyMaterialSave;
 	FL_OBJECT *MaterialProbabilitiesXYPlot;
	FL_OBJECT *Number_of_points;
	FL_OBJECT *DisplayValue;
	FL_OBJECT *InterpolateButton;
	//FL_OBJECT *MaterialTextBox;
        FL_OBJECT *material_modifiable;
        FL_OBJECT *material_count;
        void *vdata;
        long ldata;
} FD_material_sliders;
  int new_current_material;
  FL_OBJECT *obj;
  FD_material_sliders *fdui = (FD_material_sliders *) fl_calloc(1, sizeof(FD_material_sliders));
  //char which_material[20];
  double 
    oldr_ambient = *r_ambient, oldg_ambient = *g_ambient, oldb_ambient = *b_ambient,
    oldr_diffuse = *r_diffuse, oldg_diffuse = *g_diffuse, oldb_diffuse = *b_diffuse,
    oldr_specular = *r_specular, oldg_specular = *g_specular, oldb_specular = *b_specular,
    old_shinyness = *shinyness;
  float x, y; int i; char buf[64]; /* for displaying xyplot values */
  fdui->modify_material_form = fl_bgn_form(FL_NO_BOX, 610, 430);
  obj = fl_add_box(FL_UP_BOX,0,0,610,430,"");
  fdui->MaterialPropertiesTextBox = obj = fl_add_text(FL_NORMAL_TEXT,210,20,200,30,"Material Properties");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  fdui->AmbientTextBox = obj = fl_add_text(FL_NORMAL_TEXT,40,60,100,20,"ambient light");
    fl_set_object_boxtype(obj,FL_SHADOW_BOX);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->DiffuseTextBox = obj = fl_add_text(FL_NORMAL_TEXT,210,60,116,20,"diffuse reflection");
    fl_set_object_boxtype(obj,FL_SHADOW_BOX);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->SpecularTextBox = obj = fl_add_text(FL_NORMAL_TEXT,385,60,123,20,"specular reflection");
    fl_set_object_boxtype(obj,FL_SHADOW_BOX);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->ShinynessTextBox = obj = fl_add_text(FL_NORMAL_TEXT,40,228,95,20,"shinyness");
    fl_set_object_boxtype(obj,FL_SHADOW_BOX);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->AmbientRedSlider = obj = fl_add_valslider(FL_HOR_SLIDER,10,96,156,22,"Red");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_RED,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetAmbientRed,0);
  fdui->AmbientGreenSlider = obj = fl_add_valslider(FL_HOR_SLIDER,10,138,156,22,"Green");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_GREEN,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetAmbientGreen,0);
  fdui->AmbientBlueSlider = obj = fl_add_valslider(FL_HOR_SLIDER,10,179,156,22,"Blue");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_BLUE,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetAmbientBlue,0);
  fdui->DiffuseRedSlider = obj = fl_add_valslider(FL_HOR_SLIDER,187,96,156,22,"Red");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_RED,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetDiffuseRed,0);
  fdui->DiffuseGreenSlider = obj = fl_add_valslider(FL_HOR_SLIDER,187,138,156,22,"Green");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_GREEN,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetDiffuseGreen,0);
  fdui->DiffuseBlueSlider = obj = fl_add_valslider(FL_HOR_SLIDER,187,179,156,22,"Blue");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_BLUE,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetDiffuseBlue,0);
  fdui->SpecularRedSlider = obj = fl_add_valslider(FL_HOR_SLIDER,364,96,157,22,"Red");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_RED,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetSpecularRed,0);
  fdui->SpecularGreenSlider = obj = fl_add_valslider(FL_HOR_SLIDER,364,138,157,22,"Green");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_GREEN,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetSpecularGreen,0);
  fdui->SpecularBlueSlider = obj = fl_add_valslider(FL_HOR_SLIDER,364,179,157,22,"Blue");
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_BLUE,FL_LEFT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetSpecularBlue,0);
  fdui->ShinynessSlider = obj = fl_add_valslider(FL_HOR_SLIDER,10,255,156,21,"");
    fl_set_slider_bounds(obj,0.0,20.0);
    fl_set_slider_step(obj,1.0);
    fl_set_slider_return(obj, 0);
    fl_set_object_color(obj,FL_INACTIVE,FL_RIGHT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,SetShinyness,0);
  //fdui->ModifyMaterialReset = obj = fl_add_button(FL_NORMAL_BUTTON,541,20,50,20,"Reset");
  fdui->ModifyMaterialSave = obj = fl_add_button(FL_NORMAL_BUTTON,541,110,50,20,"Save");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->ModifyMaterialLoad = obj = fl_add_button(FL_NORMAL_BUTTON,541,160,50,20,"Load");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,ModifyMaterialResetCB,0);
  fdui->ModifyMaterialApply = obj = fl_add_button(FL_NORMAL_BUTTON,541,60,50,20,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  //  fl_set_object_callback(obj,ModifyMaterialApplyCB,0);
  fdui->ModifyMaterialQuit = obj = fl_add_button(FL_NORMAL_BUTTON,541,210,50,20,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
   // fl_set_object_callback(obj,ModifyMaterialQuitCB,0);
  fdui->MaterialProbabilitiesXYPlot = obj = fl_add_xyplot(FL_ACTIVE_XYPLOT,193,228,287,172,"Material Probabilities over Scalar Values");
    fl_set_object_boxtype(obj,FL_FRAME_BOX);
    fl_set_xyplot_data(obj, WX, WY, material_weight_points, "","scalar value","probability");
    fl_set_xyplot_xtics(obj, 5, 2);
    fl_set_xyplot_ytics(obj, 5, 2);
    fl_set_xyplot_ybounds(obj, 0.0, 1.0);
    fl_set_xyplot_xbounds(obj, 0.0, VP_SCALAR_MAX);
  fdui->Number_of_points = obj = fl_add_counter(FL_SIMPLE_COUNTER,510,256,69,21,"Number of points");
    fl_set_counter_precision(obj,0);
    fl_set_counter_bounds(obj,0.0,VP_SCALAR_MAX);
    fl_set_counter_value(obj,(double)material_weight_points);
    fl_set_counter_return(obj, 0);
    fl_set_object_color(obj,FL_GRAY63,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->DisplayValue = obj = fl_add_box(FL_BORDER_BOX,493,351,110,20,"");
  fdui->material_modifiable = obj = fl_add_counter(FL_SIMPLE_COUNTER,10,306,156,20,"Material for modification");
    fl_set_counter_precision(obj,0);
    fl_set_counter_bounds(obj,1.0,(double)*number_of_materials);
    fl_set_counter_value(obj,(double)*current_material);
    fl_set_counter_return(obj, 0);
    //fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_color(obj,FL_TOP_BCOL,FL_RIGHT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->material_count = obj = fl_add_counter(FL_SIMPLE_COUNTER,10,356,156,20,"Number of materials");
    fl_set_counter_precision(obj,0);
    fl_set_counter_bounds(obj,1.0,6.0);
    fl_set_counter_value(obj,(double)*number_of_materials);
    fl_set_counter_return(obj, 0);
    //fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_color(obj,FL_TOP_BCOL,FL_RIGHT_BCOL);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fl_end_form();
  fl_set_slider_value(fdui->AmbientRedSlider,oldr_ambient);
  fl_set_slider_value(fdui->AmbientGreenSlider,oldg_ambient);
  fl_set_slider_value(fdui->AmbientBlueSlider,oldb_ambient);
  fl_set_slider_value(fdui->DiffuseRedSlider,oldr_diffuse);
  fl_set_slider_value(fdui->DiffuseGreenSlider,oldg_diffuse);
  fl_set_slider_value(fdui->DiffuseBlueSlider,oldb_diffuse);
  fl_set_slider_value(fdui->SpecularRedSlider,oldr_specular);
  fl_set_slider_value(fdui->SpecularGreenSlider,oldg_specular);
  fl_set_slider_value(fdui->SpecularBlueSlider,oldb_specular);
  fl_set_slider_value(fdui->ShinynessSlider,old_shinyness);
  /* show the form */
  /*fl_show_form(fdui->modify_material_form,FL_PLACE_MOUSE,FL_TRANSIENT,"Modify Material");*/
  fl_show_form(fdui->modify_material_form,FL_PLACE_MOUSE|FL_FREE_SIZE,FL_FULLBORDER,"Modify Material");
  while (obj != fdui->ModifyMaterialQuit){
    obj = fl_do_forms();
    if (obj == fdui->Number_of_points){
      int new_material_weight_points;
      new_material_weight_points = (int)(fl_get_counter_value(obj)+0.499);
      if (new_material_weight_points != material_weight_points){
        fl_get_xyplot_data(fdui->MaterialProbabilitiesXYPlot, WX, WY, &material_weight_points);
        SetWeightPointsNumberCB(obj, material_weight_points, new_material_weight_points);
        fl_set_xyplot_data(fdui->MaterialProbabilitiesXYPlot, WX, WY, new_material_weight_points, "","scalar value","probability");
        set_material_weight_points(new_material_weight_points,*current_material);
        material_weight_points = new_material_weight_points;
      }
    }
    if (obj == fdui->material_modifiable){
      new_current_material = fl_get_counter_value(fdui->material_modifiable)+0.499; 
      if (new_current_material >  *number_of_materials) new_current_material = *number_of_materials;
      if (new_current_material != *current_material){
        *current_material = new_current_material;
        update_material_properties(
          r_ambient, g_ambient, b_ambient,
          r_diffuse, g_diffuse, b_diffuse,
          r_specular, g_specular, b_specular,
          shinyness, new_current_material, &material_weight_points);
        oldr_ambient = *r_ambient; oldg_ambient = *g_ambient; oldb_ambient = *b_ambient;
        oldr_diffuse = *r_diffuse; oldg_diffuse = *g_diffuse; oldb_diffuse = *b_diffuse;
        oldr_specular = *r_specular; oldg_specular = *g_specular; oldb_specular = *b_specular;
        fl_set_slider_value(fdui->AmbientRedSlider,oldr_ambient);
        fl_set_slider_value(fdui->AmbientGreenSlider,oldg_ambient);
        fl_set_slider_value(fdui->AmbientBlueSlider,oldb_ambient);
        fl_set_slider_value(fdui->DiffuseRedSlider,oldr_diffuse);
        fl_set_slider_value(fdui->DiffuseGreenSlider,oldg_diffuse);
        fl_set_slider_value(fdui->DiffuseBlueSlider,oldb_diffuse);
        fl_set_slider_value(fdui->SpecularRedSlider,oldr_specular);
        fl_set_slider_value(fdui->SpecularGreenSlider,oldg_specular);
        fl_set_slider_value(fdui->SpecularBlueSlider,oldb_specular);
        fl_set_slider_value(fdui->ShinynessSlider,old_shinyness);
        fl_set_xyplot_data(fdui->MaterialProbabilitiesXYPlot, WX, WY, 
               material_weight_points, "","scalar value","probability");
        fl_set_counter_value(fdui->Number_of_points,(double)material_weight_points);
      }
    } 
    if (obj == fdui->material_count){
      if (*number_of_materials != (int)(fl_get_counter_value(fdui->material_count)+0.499)){
        if ((fl_get_counter_value(fdui->material_count)+0.499) < 
            (fl_get_counter_value(fdui->material_modifiable)+0.499)){
          show_message("You can not change the total number of materials",
                     "to a value less than the current material!","");
          *number_of_materials = (int)(fl_get_counter_value(fdui->material_modifiable)+0.499);
          fl_set_counter_value(fdui->material_count,(double)*number_of_materials);
        }
        else{
          *number_of_materials = (int)(fl_get_counter_value(fdui->material_count)+0.499);
          printf("Setting new material counts to %d\n",*number_of_materials);
          fl_set_counter_bounds(fdui->material_modifiable,1.0,(double)*number_of_materials);
          update_material_numbers(*number_of_materials);
        }
      }
    }
    if (obj == fdui->MaterialProbabilitiesXYPlot){
      fl_get_xyplot(obj, &x, &y, &i);
      if (i > 0) {
        sprintf(buf,"X=%.1f  Y=%.1f",x,y);
        fl_set_object_label(fdui->DisplayValue, buf);
      }
    }
    if (obj == fdui->ModifyMaterialApply){
      FORCEDRAW = 1;
      *r_ambient = fl_get_slider_value(fdui->AmbientRedSlider); 
      *g_ambient = fl_get_slider_value(fdui->AmbientGreenSlider);
      *b_ambient = fl_get_slider_value(fdui->AmbientBlueSlider);
      *r_diffuse = fl_get_slider_value(fdui->DiffuseRedSlider); 
      *g_diffuse = fl_get_slider_value(fdui->DiffuseGreenSlider);
      *b_diffuse = fl_get_slider_value(fdui->DiffuseBlueSlider);
      *r_specular = fl_get_slider_value(fdui->SpecularRedSlider); 
      *g_specular = fl_get_slider_value(fdui->SpecularGreenSlider);
      *b_specular = fl_get_slider_value(fdui->SpecularBlueSlider);
      *shinyness = fl_get_slider_value(fdui->ShinynessSlider);
      fl_get_xyplot_data(fdui->MaterialProbabilitiesXYPlot, WX, WY, &material_weight_points);
      set_material_ramps();
      //set_material_ramp(*material_weight_points);
      modify_material_properties_applyCB(*r_ambient,*g_ambient,*b_ambient,
		       *r_diffuse,*g_diffuse,*b_diffuse,
                       *r_specular,*g_specular,*b_specular,
		       *shinyness);
    }
    if (obj == fdui->ModifyMaterialSave){
      fl_set_fselector_placement(FL_PLACE_FREE);
      //fl_set_fselector_callback(save_material_properties);
      filename = fl_show_fselector("Save Materials File", 0, "*.mat",0);
      if (filename != NULL){
        sprintf(material_filename,"%s",filename); 
        printf("Material filename is now %s\n",material_filename);
        store_materials(material_filename);
      }
      else{
        show_message("No file name selected \n","","");
      }
    }
    if (obj == fdui->ModifyMaterialLoad){
      fl_set_fselector_placement(FL_PLACE_FREE);
      //fl_set_fselector_callback(load_material_properties);
      filename = fl_show_fselector("Load Materials File", 0, "*.mat",0);
      if (filename != NULL){
        sprintf(material_filename,"%s",filename); 
        printf("Material filename is now %s\n",material_filename);
        load_materials(material_filename);
        new_current_material = fl_get_counter_value(fdui->material_modifiable)+0.499; 
        update_material_properties(
          r_ambient, g_ambient, b_ambient,
          r_diffuse, g_diffuse, b_diffuse,
          r_specular, g_specular, b_specular,
          shinyness, new_current_material, &material_weight_points);
        *number_of_materials = get_number_of_materials();
        printf("number_of_materials is now %d\n",*number_of_materials);
        if (*number_of_materials != (int)(fl_get_counter_value(fdui->material_count)+0.499)){
          fl_set_counter_value(fdui->material_count,(double)*number_of_materials);
          fl_set_counter_bounds(fdui->material_modifiable,1.0,(double)*number_of_materials);
          update_material_numbers(*number_of_materials);
        }
        oldr_ambient = *r_ambient; oldg_ambient = *g_ambient; oldb_ambient = *b_ambient;
        oldr_diffuse = *r_diffuse; oldg_diffuse = *g_diffuse; oldb_diffuse = *b_diffuse;
        oldr_specular = *r_specular; oldg_specular = *g_specular; oldb_specular = *b_specular;
        fl_set_slider_value(fdui->AmbientRedSlider,oldr_ambient);
        fl_set_slider_value(fdui->AmbientGreenSlider,oldg_ambient);
        fl_set_slider_value(fdui->AmbientBlueSlider,oldb_ambient);
        fl_set_slider_value(fdui->DiffuseRedSlider,oldr_diffuse);
        fl_set_slider_value(fdui->DiffuseGreenSlider,oldg_diffuse);
        fl_set_slider_value(fdui->DiffuseBlueSlider,oldb_diffuse);
        fl_set_slider_value(fdui->SpecularRedSlider,oldr_specular);
        fl_set_slider_value(fdui->SpecularGreenSlider,oldg_specular);
        fl_set_slider_value(fdui->SpecularBlueSlider,oldb_specular);
        fl_set_slider_value(fdui->ShinynessSlider,old_shinyness);
        fl_set_xyplot_data(fdui->MaterialProbabilitiesXYPlot, WX, WY, 
               material_weight_points, "","scalar value","probability");
        fl_set_counter_value(fdui->Number_of_points,(double)material_weight_points);
        FORCEDRAW = 1;
        modify_material_properties_applyCB(*r_ambient,*g_ambient,*b_ambient,
		       *r_diffuse,*g_diffuse,*b_diffuse,
                       *r_specular,*g_specular,*b_specular,
		       *shinyness);
      }
      else{
        show_message("No file name selected! \n","","");
      }
    }
    //if (obj == fdui->ModifyMaterialReset){
    //  fl_set_slider_value(fdui->AmbientRedSlider,oldr_ambient); *r_ambient = oldr_ambient;
    //  fl_set_slider_value(fdui->AmbientGreenSlider,oldg_ambient); *g_ambient = oldg_ambient;
    //  fl_set_slider_value(fdui->AmbientBlueSlider,oldb_ambient); *b_ambient = oldb_ambient;
    //  fl_set_slider_value(fdui->DiffuseRedSlider,oldr_diffuse); *r_diffuse = oldr_diffuse;
    //  fl_set_slider_value(fdui->DiffuseGreenSlider,oldg_diffuse); *g_diffuse = oldg_diffuse;
    //  fl_set_slider_value(fdui->DiffuseBlueSlider,oldb_diffuse); *b_diffuse = oldb_diffuse;
    //  fl_set_slider_value(fdui->SpecularRedSlider,oldr_specular); *r_specular = oldr_specular;
    //  fl_set_slider_value(fdui->SpecularGreenSlider,oldg_specular); *g_specular = oldg_specular;
    //  fl_set_slider_value(fdui->SpecularBlueSlider,oldb_specular); *b_specular = oldb_specular;
    //  fl_set_slider_value(fdui->ShinynessSlider,old_shinyness); *shinyness = old_shinyness;
    //}
  }
  fl_hide_form(fdui->modify_material_form);
  fl_free_form(fdui->modify_material_form);
}

void depth_cueing_slider(double *front_factor, double *fog_density, char enabled)
/********************************************************************************/
/* Activates the depth cueing pulldown selection from the menu modify   	*/
/* For changing the depth cueing values front_factor and fog_density where the  */
/* fraction of light transmitted after travelling a distance d through a fog is */
/* proportional to the front_factor times exp(-fog_density*d)                   */
/********************************************************************************/
{
  typedef struct {
        FL_FORM *depth_cueing_slider;
        FL_OBJECT *front_factor;
        FL_OBJECT *DepthCueApply;
        FL_OBJECT *DepthCueReset;
        FL_OBJECT *DepthCueQuit;
        FL_OBJECT *fog_density;
        FL_OBJECT *DepthCueingTextBox;
        FL_OBJECT *InformationTextBox1;
        FL_OBJECT *InformationTextBox2;
        void *vdata;
        long ldata;
  } FD_depth_cueing_slider;
  FL_OBJECT *obj;
  FD_depth_cueing_slider *fdui = (FD_depth_cueing_slider *)  fl_calloc(1, sizeof(FD_depth_cueing_slider));
  double old_front_factor = *front_factor, old_fog_density = *fog_density;
  fdui->depth_cueing_slider = fl_bgn_form(FL_NO_BOX, 300, 320);
  obj = fl_add_box(FL_UP_BOX,0,0,300,320,"");
  fdui->front_factor = obj = fl_add_valslider(FL_HOR_SLIDER,30,130,230,30,"Front factor");
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.0,10.0);
    fl_set_slider_step(obj,0.05);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_value(fdui->front_factor,*front_factor); 
  fdui->fog_density = obj = fl_add_valslider(FL_HOR_SLIDER,30,190,230,30,"Fog density");
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.0,20.0);
    fl_set_slider_step(obj,0.2);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_value(fdui->fog_density,*fog_density); 
  fdui->DepthCueApply = obj = fl_add_button(FL_NORMAL_BUTTON,30,260,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->DepthCueReset = obj = fl_add_button(FL_NORMAL_BUTTON,110,260,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->DepthCueQuit = obj = fl_add_button(FL_NORMAL_BUTTON,190,260,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->DepthCueingTextBox = obj = fl_add_text(FL_NORMAL_TEXT,70,20,170,30,"Depth Cueing");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  if (enabled){
    fdui->InformationTextBox1 = obj = fl_add_text(FL_NORMAL_TEXT,30,60,230,20,"currently enabled");
    fdui->InformationTextBox2 = obj = fl_add_text(FL_NORMAL_TEXT,30,90,230,20,"to disable use switch menu ");
  }
  else{
    fdui->InformationTextBox1 = obj = fl_add_text(FL_NORMAL_TEXT,30,60,230,20,"currently not enabled");
    fdui->InformationTextBox2 = obj = fl_add_text(FL_NORMAL_TEXT,30,90,230,20,"to enable use switch menu ");
  }
  fl_end_form();
  fl_show_form(fdui->depth_cueing_slider,FL_PLACE_MOUSE,FL_TRANSIENT,"Modify Depth Cueing");
  while (obj != fdui->DepthCueQuit){
    obj = fl_do_forms();
    if (obj == fdui->DepthCueApply){
      *front_factor = fl_get_slider_value(fdui->front_factor); 
      *fog_density = fl_get_slider_value(fdui->fog_density);
      setDepthCueingValues(*front_factor,*fog_density);
    } 
    if (obj == fdui->DepthCueReset){
      fl_set_slider_value(fdui->front_factor,old_front_factor); *front_factor = old_front_factor;
      fl_set_slider_value(fdui->fog_density,old_fog_density); *fog_density = old_fog_density;
      setDepthCueingValues(*front_factor,*fog_density);
    }
  }
  fl_hide_form(fdui->depth_cueing_slider);
  fl_free_form(fdui->depth_cueing_slider);
}

void maximum_opacity_slider(double *maximum_opacity)
/********************************************************************************/
/* Activates the maximum opacity pulldown selection from the menu modify	*/
/* maximum opacity sets the opacity level for total occlusion of anything 	*/
/* along a ray behind a voxel if this accumulated opacity is reached 		*/
/********************************************************************************/
{
  typedef struct {
        FL_FORM *maximum_opacity;
        FL_OBJECT *MaxOpacitySlider;
        FL_OBJECT *MaxOpacityApply;
        FL_OBJECT *MaxOpacityReset;
        FL_OBJECT *MaxOpacityQuit;
        FL_OBJECT *MaxOpacityTextBox;
        FL_OBJECT *InformationTextBox1;
        FL_OBJECT *InformationTextBox2;
        void *vdata;
        long ldata;
  } FD_maximum_opacity_slider;

  FL_OBJECT *obj;
  FD_maximum_opacity_slider *fdui = 
   (FD_maximum_opacity_slider *)  fl_calloc(1, sizeof(FD_maximum_opacity_slider));
  double old_maximum_opacity = *maximum_opacity;
  fdui->maximum_opacity = fl_bgn_form(FL_NO_BOX, 310, 260);
  obj = fl_add_box(FL_UP_BOX,0,0,310,260,"");
  fdui->MaxOpacitySlider = obj = fl_add_valslider(FL_HOR_SLIDER,30,130,250,30,"Maximum Opacity");
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.0,1.0);
    fl_set_slider_step(obj,0.01);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_value(fdui->MaxOpacitySlider,*maximum_opacity); 
  fdui->MaxOpacityApply = obj = fl_add_button(FL_NORMAL_BUTTON,30,200,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->MaxOpacityReset = obj = fl_add_button(FL_NORMAL_BUTTON,120,200,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->MaxOpacityQuit = obj = fl_add_button(FL_NORMAL_BUTTON,210,200,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->MaxOpacityTextBox = obj = fl_add_text(FL_NORMAL_TEXT,60,20,230,30,"Maximum Opacity");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  fdui->InformationTextBox1 = obj = fl_add_text(FL_NORMAL_TEXT,30,60,260,20,
         "Early termination of ray casting if accumulated");
  fdui->InformationTextBox2 = obj = fl_add_text(FL_NORMAL_TEXT,30,90,260,20,
         "opacity along a ray exceeds maximum opacity");
  fl_end_form();
  fl_show_form(fdui->maximum_opacity,FL_PLACE_MOUSE,FL_TRANSIENT,"Modify Maximum Opacity");
  while (obj != fdui->MaxOpacityQuit){
    obj = fl_do_forms();
    if (obj == fdui->MaxOpacityApply){
      *maximum_opacity = fl_get_slider_value(fdui->MaxOpacitySlider);
      setMaximumOpacity(*maximum_opacity); 
      printf("Apply button pushed");
      multi_rotate_image(1,0.0,1);
    } 
    if (obj == fdui->MaxOpacityReset){
      fl_set_slider_value(fdui->MaxOpacitySlider,(float)old_maximum_opacity); 
      setMaximumOpacity(old_maximum_opacity); 
      *maximum_opacity = old_maximum_opacity;
      multi_rotate_image(1,0.0,1);
    }
  }
  fl_hide_form(fdui->maximum_opacity);
  fl_free_form(fdui->maximum_opacity);
}

void minimum_opacity_slider(double *minimum_opacity)
/********************************************************************************/
/* Activates the minimum opacity pulldown selection from the menu modify	*/
/* all voxels having an opacity below the minimum opacity are considered        */
/* transparent and may be skipped						*/
/********************************************************************************/
{
  typedef struct {
        FL_FORM *minimum_opacity_form;
        FL_OBJECT *MinOpacitySlider;
        FL_OBJECT *MinOpacityApply;
        FL_OBJECT *MinOpacityReset;
        FL_OBJECT *MinOpacityQuit;
        FL_OBJECT *MinOpacityTextBox;
        FL_OBJECT *InformationTextBox1;
        FL_OBJECT *InformationTextBox2;
        void *vdata;
        long ldata;
  } FD_minimum_opacity_slider;

  FL_OBJECT *obj;
  FD_minimum_opacity_slider *fdui = 
   (FD_minimum_opacity_slider *)  fl_calloc(1, sizeof(FD_minimum_opacity_slider));
  double old_minimum_opacity = *minimum_opacity;
  fdui->minimum_opacity_form = fl_bgn_form(FL_NO_BOX, 360, 260);
  obj = fl_add_box(FL_UP_BOX,0,0,360,260,"");
  fdui->MinOpacitySlider = obj = fl_add_valslider(FL_HOR_SLIDER,30,130,300,30,"Minimum Opacity");
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.0,1.0);
    fl_set_slider_step(obj,0.01);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_value(fdui->MinOpacitySlider,*minimum_opacity); 
  fdui->MinOpacityApply = obj = fl_add_button(FL_NORMAL_BUTTON,30,200,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->MinOpacityReset = obj = fl_add_button(FL_NORMAL_BUTTON,145,200,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->MinOpacityQuit = obj = fl_add_button(FL_NORMAL_BUTTON,260,200,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->MinOpacityTextBox = obj = fl_add_text(FL_NORMAL_TEXT,90,20,230,30,"Minimum Opacity");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  fdui->InformationTextBox1 = obj = fl_add_text(FL_NORMAL_TEXT,25,60,320,20,
         "If the opacity of a voxel is less or equal to the minimum");
  fdui->InformationTextBox2 = obj = fl_add_text(FL_NORMAL_TEXT,25,90,320,20,
         "opacity, then the voxel is considered totally transparent");
  fl_end_form();
  fl_show_form(fdui->minimum_opacity_form,FL_PLACE_MOUSE,FL_TRANSIENT,"Modify Minimum Opacity");
  while (obj != fdui->MinOpacityQuit){
    obj = fl_do_forms();
    if (obj == fdui->MinOpacityApply){
      *minimum_opacity = fl_get_slider_value(fdui->MinOpacitySlider); 
      setMinimumOpacity(*minimum_opacity);
      new_classification();
      //multi_rotate_image(1,0.0,1);
    } 
    if (obj == fdui->MinOpacityReset){
      fl_set_slider_value(fdui->MinOpacitySlider,(float)old_minimum_opacity); 
      *minimum_opacity = old_minimum_opacity;
      setMinimumOpacity(*minimum_opacity);
      new_classification();
      //multi_rotate_image(1,0.0,1);
    }
  }
  fl_hide_form(fdui->minimum_opacity_form);
  fl_free_form(fdui->minimum_opacity_form);
}

void show_message(const char *s1, const char *s2, const char *s3)
/************************************************************************/
/*  displays a message							*/
/************************************************************************/
{
  fl_show_message(s1,s2,s3);
}

//int show_question(const char *s1, const char *s2, const char *s3)
int show_question(const char *s1)
/************************************************************************/
/*  forms a yes-no-choice question					*/
/************************************************************************/
{
  int return_id;
  //return_id = fl_show_question(s1,s2,s3);
  return_id = fl_show_question(s1,1);
  return return_id;
}

void ask_server_name()
/************************************************************************/
/*  requires a string from a input form					*/
/************************************************************************/
{
  typedef struct {
        FL_FORM *form;
        FL_OBJECT *Question_Server_Name;
        void *vdata;
        long ldata;
  } FD_ask_for_server;

  FL_OBJECT *obj;
  FD_ask_for_server *fdui = (FD_ask_for_server *) fl_calloc(1, sizeof(FD_ask_for_server));
  fdui->form = fl_bgn_form(FL_NO_BOX, 530, 120);
    obj = fl_add_box(FL_UP_BOX,0,0,530,120,"");
    fdui->Question_Server_Name = obj = fl_add_input(FL_NORMAL_INPUT,190,40,300,40,"Server Name");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
    fl_set_object_lstyle(obj,FL_TIMESBOLD_STYLE);
    fl_set_input(obj, server_name);
  fl_end_form();

  fl_show_form(fdui->form,FL_PLACE_MOUSE,FL_TRANSIENT,"Please enter server name");
  obj = fl_do_forms();
  sprintf(server_name,fl_get_input(obj));
  printf("server name is %s\n", server_name);
  fl_hide_form(fdui->form);
  fl_free_form(fdui->form);
}

void ask_port_address()
/************************************************************************/
/*  requires a string from a input form					*/
/************************************************************************/
{
  typedef struct {
        FL_FORM *form;
        FL_OBJECT *Question_Port_Address;
        void *vdata;
        long ldata;
  } FD_ask_for_port;
  char base_port_ascii[10];

  FL_OBJECT *obj;
  sprintf(base_port_ascii,"%d",base_port);
  printf("old port address is %s\n", base_port_ascii);
  FD_ask_for_port *fdui = (FD_ask_for_port *) fl_calloc(1, sizeof(FD_ask_for_port));
  fdui->form = fl_bgn_form(FL_NO_BOX, 530, 120);
    obj = fl_add_box(FL_UP_BOX,0,0,530,120,"");
    fdui->Question_Port_Address = obj = fl_add_input(FL_NORMAL_INPUT,190,40,300,40,"Port Address");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
    fl_set_object_lstyle(obj,FL_TIMESBOLD_STYLE);
    fl_set_input(obj, base_port_ascii);
  fl_end_form();

  fl_show_form(fdui->form,FL_PLACE_MOUSE,FL_TRANSIENT,"Please enter new port address");
  obj = fl_do_forms();
  sprintf(base_port_ascii,fl_get_input(obj));
  printf("new port address is %s\n", base_port_ascii);
  base_port = (short int) atoi(base_port_ascii);
  fl_hide_form(fdui->form);
  fl_free_form(fdui->form);
}

void set_dimensions(int *data_header, int *data_xlen, int *data_ylen, 
                    int *data_zlen, unsigned char *apply_button_pushed)
/********************************************************************************/
/* input of the volume sizes in the three dimensions and the header size	*/
/* to be skipped is required if a new unclassified volume is loaded		*/
/********************************************************************************/
{

  typedef struct {
        FL_FORM *set_dimensions_form;
        FL_OBJECT *SetDimensionsTextBox;
        FL_OBJECT *InformationTextBox1;
        FL_OBJECT *InformationTextBox2;
	FL_OBJECT *HeaderSizeInput;
        FL_OBJECT *SizeXInput;
        FL_OBJECT *SizeYInput;
        FL_OBJECT *SizeZInput;
        FL_OBJECT *SetDimApply;
        FL_OBJECT *SetDimReset;
        FL_OBJECT *SetDimQuit;
        void *vdata;
        long ldata;
  } FD_set_dimensions;

  int old_data_header, old_data_xlen, old_data_ylen, old_data_zlen;
  old_data_header = *data_header;
  old_data_xlen = *data_xlen;
  old_data_ylen = *data_ylen;
  old_data_zlen = *data_zlen;
  char bufhead[8], bufx[8], bufy[8], bufz[8];
  sprintf(bufhead,"%d",old_data_header);
  sprintf(bufx,"%d",old_data_xlen);
  sprintf(bufy,"%d",old_data_ylen);
  sprintf(bufz,"%d",old_data_zlen);
 

  FL_OBJECT *obj;
  FD_set_dimensions *fdui = (FD_set_dimensions *) 
     fl_calloc(1, sizeof(FD_set_dimensions));
  fdui->set_dimensions_form = fl_bgn_form(FL_NO_BOX, 300, 410);
    obj = fl_add_box(FL_UP_BOX,0,0,300,410,"");
  fdui->SetDimensionsTextBox = obj = fl_add_text(FL_NORMAL_TEXT,70,20,170,30,"Set Dimensions");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  fdui->InformationTextBox1 = obj = fl_add_text(FL_NORMAL_TEXT,30,60,230,20,"input the header size");
    fl_set_object_lalign(obj,FL_ALIGN_CENTER);
  fdui->HeaderSizeInput = obj = fl_add_input(FL_NORMAL_INPUT,140,90,50,30,"header size");
    fl_set_input(obj,bufhead);
  fdui->InformationTextBox2 = obj = fl_add_text(FL_NORMAL_TEXT,30,140,230,20,"input of the volume  dimensions");
    fl_set_object_lalign(obj,FL_ALIGN_CENTER);
  fdui->SizeXInput = obj = fl_add_input(FL_INT_INPUT,140,180,40,30,"size x");
    fl_set_input(obj,bufx);
  fdui->SizeYInput = obj = fl_add_input(FL_INT_INPUT,140,230,40,30,"size y");
    fl_set_input(obj,bufy);
  fdui->SizeZInput = obj = fl_add_input(FL_INT_INPUT,140,280,40,30,"size z");
    fl_set_input(obj,bufz);
  fdui->SetDimApply = obj = fl_add_button(FL_NORMAL_BUTTON,30,340,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->SetDimReset = obj = fl_add_button(FL_NORMAL_BUTTON,110,340,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->SetDimQuit = obj = fl_add_button(FL_NORMAL_BUTTON,190,340,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fl_end_form();
  obj = fdui->InformationTextBox1;
  fl_show_form(fdui->set_dimensions_form,FL_PLACE_MOUSE,FL_TRANSIENT,"Set dimensions");
  while (obj != fdui->SetDimQuit){
    obj = fl_do_forms();
    if (obj == fdui->SetDimApply){
      *data_header = atoi(fl_get_input(fdui->HeaderSizeInput));
      *data_xlen = atoi(fl_get_input(fdui->SizeXInput)); 
      *data_ylen = atoi(fl_get_input(fdui->SizeYInput)); 
      *data_zlen = atoi(fl_get_input(fdui->SizeZInput)); 
      *apply_button_pushed = 1;
      obj = fdui->SetDimQuit;
    } 
    if (obj == fdui->SetDimReset){
       *data_header = old_data_header;
       *data_xlen = old_data_xlen;
       *data_ylen = old_data_ylen;
       *data_zlen = old_data_zlen;
       sprintf(bufhead,"%d",old_data_header);
       sprintf(bufx,"%d",old_data_xlen);
       sprintf(bufy,"%d",old_data_ylen);
       sprintf(bufz,"%d",old_data_zlen);
       fl_set_input(fdui->HeaderSizeInput,bufhead);
       fl_set_input(fdui->SizeXInput,bufx);
       fl_set_input(fdui->SizeYInput,bufy);
       fl_set_input(fdui->SizeZInput,bufz);
    }
  }
  fl_hide_form(fdui->set_dimensions_form);
  fl_free_form(fdui->set_dimensions_form);
}

void scale_voxels(double *scale_x, double *scale_y, double *scale_z)
/************************************************************************/
/* Activates the scaling pulldown selection from the menu modify	*/
/* for scaling the relative distances between the voxels		*/
/************************************************************************/
{
typedef struct {
        FL_FORM *scale_form;
        FL_OBJECT *ScaleX;
        FL_OBJECT *ScaleY;
        FL_OBJECT *ScaleZ;
        FL_OBJECT *ScaleApply;
        FL_OBJECT *ScaleReset;
        FL_OBJECT *ScaleQuit;
        FL_OBJECT *ScalingTextBox;
        FL_OBJECT *InformationTextBox1;
        void *vdata;
        long ldata;
} FD_scale_sliders;

  FL_OBJECT *obj;
  FD_scale_sliders *fdui = (FD_scale_sliders *) 
     fl_calloc(1, sizeof(FD_scale_sliders));


  fdui->scale_form = fl_bgn_form(FL_NO_BOX, 300, 330);
  obj = fl_add_box(FL_UP_BOX,0,0,300,330,"");
  fdui->ScaleX = obj = fl_add_valslider(FL_HOR_SLIDER,30,90,230,30,"Scale x");
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.1,20.0);
    fl_set_slider_step(obj,0.1);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_value(fdui->ScaleX,*scale_x); 
  fdui->ScaleY = obj = fl_add_valslider(FL_HOR_SLIDER,30,150,230,30,"Scale y");
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.1,20.0);
    fl_set_slider_step(obj,0.1);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_value(fdui->ScaleY,*scale_y); 
  fdui->ScaleZ = obj = fl_add_valslider(FL_HOR_SLIDER,30,210,230,30,"Scale z");
    fl_set_slider_return(obj, 0);
    fl_set_slider_bounds(obj,0.1,20.0);
    fl_set_slider_step(obj,0.1);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_slider_value(fdui->ScaleZ,*scale_z); 
  fdui->ScaleApply = obj = fl_add_button(FL_NORMAL_BUTTON,30,280,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->ScaleReset = obj = fl_add_button(FL_NORMAL_BUTTON,110,280,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->ScaleQuit = obj = fl_add_button(FL_NORMAL_BUTTON,190,280,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->ScalingTextBox = obj = fl_add_text(FL_NORMAL_TEXT,70,20,170,30,"Scaling");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  fdui->InformationTextBox1 = obj = fl_add_text(FL_NORMAL_TEXT,30,60,230,20,"sets the relative  distances between voxels");
  fl_end_form();
  fl_show_form(fdui->scale_form,FL_PLACE_MOUSE,FL_TRANSIENT,"Scale Values");
  while (obj != fdui->ScaleQuit){
    obj = fl_do_forms();
    if (obj == fdui->ScaleApply){
      *scale_x = fl_get_slider_value(fdui->ScaleX); 
      *scale_y = fl_get_slider_value(fdui->ScaleY); 
      *scale_z = fl_get_slider_value(fdui->ScaleZ);
      multi_rotate_image(1,0.0,1); 
    } 
    if (obj == fdui->ScaleReset){
      fl_set_slider_value(fdui->ScaleX,1.0); 
      fl_set_slider_value(fdui->ScaleY,1.0); 
      fl_set_slider_value(fdui->ScaleZ,1.0);
      *scale_x = 1.0;
      *scale_y = 1.0;
      *scale_z = 1.0; 
    }
  }
  fl_hide_form(fdui->scale_form);
  fl_free_form(fdui->scale_form);
}

void material_numbers_sliders(int *material_count, int *current_material)
/********************************************************************************/
/* Activates the material numbers pulldown selection from the menu modify	*/
/* for modifying material_count and current_material where			*/
/* current_material is the material that can be modified in the pulldown menu	*/
/* material properties and material_count is the total number of materials used */
/********************************************************************************/
{
typedef struct {
        FL_FORM *material_numbers_form;
        FL_OBJECT *MaterialNumberTextBox;
        FL_OBJECT *InformationTextBox1;
        FL_OBJECT *InformationTextBox2;
        FL_OBJECT *InformationTextBox3;
        FL_OBJECT *material_modifiable;
        FL_OBJECT *material_count;
        FL_OBJECT *MaterialNumberApply;
        FL_OBJECT *MaterialNumberReset;
        FL_OBJECT *MaterialNumberQuit;
        void *vdata;
        long ldata;
} FD_material_numbers_sliders;
  FL_OBJECT *obj;
  FD_material_numbers_sliders *fdui = (FD_material_numbers_sliders *) 
     fl_calloc(1, sizeof(FD_material_numbers_sliders));
  int old_material_count = *material_count, old_current_material = *current_material;
  fdui->material_numbers_form = fl_bgn_form(FL_NO_BOX, 300, 310);
  obj = fl_add_box(FL_UP_BOX,0,0,300,310,"");
  fdui->MaterialNumberTextBox = obj = fl_add_text(FL_NORMAL_TEXT,70,10,170,30,"Material Number");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  fdui->material_modifiable = obj = fl_add_counter(FL_SIMPLE_COUNTER,30,100,230,30,"Material for modification");
    fl_set_counter_precision(obj,0);
    fl_set_counter_bounds(obj,1.0,(double)*material_count);
    fl_set_counter_value(obj,(double)*current_material);
    fl_set_counter_return(obj, 0);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->material_count = obj = fl_add_counter(FL_SIMPLE_COUNTER,30,190,230,30,"Number of materials");
    fl_set_counter_precision(obj,0);
    fl_set_counter_bounds(obj,1.0,6.0);
    fl_set_counter_value(obj,(double)*material_count);
    fl_set_counter_return(obj, 0);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->MaterialNumberApply = obj = fl_add_button(FL_NORMAL_BUTTON,30,260,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->MaterialNumberReset = obj = fl_add_button(FL_NORMAL_BUTTON,110,260,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->MaterialNumberQuit = obj = fl_add_button(FL_NORMAL_BUTTON,190,260,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->InformationTextBox1 = obj = fl_add_text(FL_NORMAL_TEXT,30,50,250,20,"current material number for modification");
  fdui->InformationTextBox2 = obj = fl_add_text(FL_NORMAL_TEXT,30,160,230,20,"total number of materials (maximum 6)");
  fdui->InformationTextBox3 = obj = fl_add_text(FL_NORMAL_TEXT,30,70,250,20,"in submenu material properties ...");
  fl_end_form();
  fl_show_form(fdui->material_numbers_form,FL_PLACE_MOUSE,FL_TRANSIENT,"Material Numbers ...");
  while (obj != fdui->MaterialNumberQuit){
    obj = fl_do_forms();
    if (obj == fdui->MaterialNumberApply){
      *current_material = fl_get_counter_value(fdui->material_modifiable)+0.499; 
      *material_count = fl_get_counter_value(fdui->material_count)+0.499;
      printf("current_material and material_count is: %d %d\n",*current_material,*material_count);
      fl_set_counter_bounds(fdui->material_modifiable,1.0,(double)*material_count);
      if (*current_material >  *material_count) *current_material = *material_count;
    } 
    if (obj == fdui->MaterialNumberReset){
      fl_set_counter_value(fdui->material_modifiable,(double)old_current_material); 
      *current_material = old_current_material;
      fl_set_counter_value(fdui->material_count,(double)old_material_count);
      *material_count = old_material_count; 
      fl_set_counter_bounds(fdui->material_modifiable,1.0,(double)*material_count);
    }
  }
  fl_hide_form(fdui->material_numbers_form);
  fl_free_form(fdui->material_numbers_form);
}

void light_numbers_menu(int *current_light, int *light_enabled)
/********************************************************************************/
/* Activates the light numbers pulldown selection from the menu modify	        */
/* for en-/disabling directional light sources and for selecting the 		*/
/* light source that may be modified in the pulldown menu light properties	*/
/********************************************************************************/
{
typedef struct {
	FL_FORM *light_number_form;
	FL_OBJECT *light_modifiable;
	FL_OBJECT *LightNumberApply;
	FL_OBJECT *LightNumberReset;
	FL_OBJECT *LightNumberQuit;
	FL_OBJECT *LightNumberTextBox;
	FL_OBJECT *Light1_Button;
	FL_OBJECT *Light2_Button;
	FL_OBJECT *Light3_Button;
	FL_OBJECT *Light4_Button;
	FL_OBJECT *Light5_Button;
	FL_OBJECT *Light6_Button;
	FL_OBJECT *InformationTextBox1;
	FL_OBJECT *InformationTextBox2;
	FL_OBJECT *InformationTextBox3;
	void *vdata;
	long ldata;
} FD_light_number_struct;

  int old_current_light = *current_light;
  int orig_light_enabled[6];
  FL_OBJECT *obj;
  FD_light_number_struct *fdui = (FD_light_number_struct *) fl_calloc(1, 
    sizeof(FD_light_number_struct));
  fdui->light_number_form = fl_bgn_form(FL_NO_BOX, 300, 420);
  obj = fl_add_box(FL_UP_BOX,0,0,300,420,"");
  fdui->light_modifiable = obj = fl_add_counter(FL_SIMPLE_COUNTER,30,100,230,30,"Light for modification");
    fl_set_counter_precision(obj,0);
    fl_set_counter_bounds(obj,1.0,6.0);
    fl_set_counter_value(obj,(double)*current_light);
    fl_set_counter_return(obj, 0);
    fl_set_object_color(obj,FL_BLUE,FL_RED);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->LightNumberApply = obj = fl_add_button(FL_NORMAL_BUTTON,30,360,70,30,"Apply");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->LightNumberReset = obj = fl_add_button(FL_NORMAL_BUTTON,110,360,70,30,"Reset");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->LightNumberQuit = obj = fl_add_button(FL_NORMAL_BUTTON,190,360,70,30,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fdui->LightNumberTextBox = obj = fl_add_text(FL_NORMAL_TEXT,30,10,250,30,"Enable and Select Lights");
    fl_set_object_lsize(obj,FL_LARGE_SIZE);
  fdui->Light1_Button = obj = fl_add_lightbutton(FL_PUSH_BUTTON,50,200,80,30,"Light 1");
    fl_set_button(obj, *light_enabled);
  fdui->Light2_Button = obj = fl_add_lightbutton(FL_PUSH_BUTTON,50,250,80,30,"Light 2");
    fl_set_button(obj, *(light_enabled+1));
  fdui->Light3_Button = obj = fl_add_lightbutton(FL_PUSH_BUTTON,50,300,80,30,"Light 3");
    fl_set_button(obj, *(light_enabled+2));
  fdui->Light4_Button = obj = fl_add_lightbutton(FL_PUSH_BUTTON,160,200,80,30,"Light 4");
    fl_set_button(obj, *(light_enabled+3));
  fdui->Light5_Button = obj = fl_add_lightbutton(FL_PUSH_BUTTON,160,250,80,30,"Light 5");
    fl_set_button(obj, *(light_enabled+4));
  fdui->Light6_Button = obj = fl_add_lightbutton(FL_PUSH_BUTTON,160,300,80,30,"Light 6");
    fl_set_button(obj, *(light_enabled+5));
  fdui->InformationTextBox1 = obj = fl_add_text(FL_NORMAL_TEXT,30,50,250,20,"current light number for modification");
  fdui->InformationTextBox2 = obj = fl_add_text(FL_NORMAL_TEXT,30,70,250,20,"in submenu light properties ...");
  fdui->InformationTextBox3 = obj = fl_add_text(FL_NORMAL_TEXT,30,170,250,20,"enable / disable  lights");
  fl_end_form();
  *orig_light_enabled = *light_enabled;
  *(orig_light_enabled+1) = *(light_enabled+1);
  *(orig_light_enabled+2) = *(light_enabled+2);
  *(orig_light_enabled+3) = *(light_enabled+3);
  *(orig_light_enabled+4) = *(light_enabled+4);
  *(orig_light_enabled+5) = *(light_enabled+5);
  fl_show_form(fdui->light_number_form,FL_PLACE_MOUSE,FL_TRANSIENT,"Light Numbers ...");
  while (obj != fdui->LightNumberQuit){
    obj = fl_do_forms();
    if (obj == fdui->LightNumberApply){
      *current_light = fl_get_counter_value(fdui->light_modifiable)+0.499;
      *(light_enabled) = fl_get_button(fdui->Light1_Button); 
      *(light_enabled+1) = fl_get_button(fdui->Light2_Button); 
      *(light_enabled+2) = fl_get_button(fdui->Light3_Button); 
      *(light_enabled+3) = fl_get_button(fdui->Light4_Button); 
      *(light_enabled+4) = fl_get_button(fdui->Light5_Button); 
      *(light_enabled+5) = fl_get_button(fdui->Light6_Button); 
      setLights(light_enabled);
      multi_rotate_image(1,0.0,1); 
      //printf("current_light is: %d \n",*current_light);
    } 
    if (obj == fdui->LightNumberReset){
      fl_set_counter_value(fdui->light_modifiable,(double)old_current_light); 
      *current_light = old_current_light;
      *light_enabled = *orig_light_enabled;
      *(light_enabled+1) = *(orig_light_enabled+1);
      *(light_enabled+2) = *(orig_light_enabled+2);
      *(light_enabled+3) = *(orig_light_enabled+3);
      *(light_enabled+4) = *(orig_light_enabled+4);
      *(light_enabled+5) = *(orig_light_enabled+5);
      fl_set_button(fdui->Light1_Button, *light_enabled);
      fl_set_button(fdui->Light2_Button, *(light_enabled+1));
      fl_set_button(fdui->Light3_Button, *(light_enabled+2));
      fl_set_button(fdui->Light4_Button, *(light_enabled+3));
      fl_set_button(fdui->Light5_Button, *(light_enabled+4));
      fl_set_button(fdui->Light6_Button, *(light_enabled+5));
      setLights(light_enabled);
      multi_rotate_image(1,0.0,1); 
    }
  }
  fl_hide_form(fdui->light_number_form);
  fl_free_form(fdui->light_number_form);
}


FD_histogram * create_histogram_object(){
/**********************************************/
/* create the histogram form stuff	      */
/**********************************************/
  FD_histogram *fdui;
  FL_OBJECT *obj;
  float last_received_minimum, last_received_maximum;
  char string_for_minimum[50], string_for_maximum[50];
  
  last_received_minimum = 0;
  last_received_maximum = 0;
  fdui = (FD_histogram *) fl_calloc(1, sizeof(FD_histogram));
  fdui->histogram_form = fl_bgn_form(FL_NO_BOX, 560, 360);
    obj = fl_add_box(FL_UP_BOX,0,0,560,360,"");
    fdui->histogram_xyplot = obj = fl_add_xyplot(FL_IMPULSE_XYPLOT,40,50,460,220,"Histogram");
    fl_set_object_boxtype(obj,FL_UP_BOX);
  fdui->MinimumHint = obj = fl_add_text(FL_NORMAL_TEXT,40,310,230,20,string_for_minimum);
    fl_set_object_boxtype(obj,FL_UP_BOX);
  fdui->MaximumHint = obj = fl_add_text(FL_NORMAL_TEXT,270,310,230,20,string_for_maximum);
    fl_set_object_boxtype(obj,FL_UP_BOX);
    //fdui->histogram_quit = obj = fl_add_button(FL_NORMAL_BUTTON,430,310,70,30,"Quit");
  fl_end_form();
  return fdui;
}

void remove_histogram(){
/**********************************************/
/* removes the histogram     		      */
/**********************************************/
  fl_hide_form(histogram_object->histogram_form);
  fl_free_form(histogram_object->histogram_form);
}

void draw_histogram(unsigned int* histogram_array, int histogram_min, int histogram_max){
/**********************************************/
/* draws the histogram     		      */
/**********************************************/
  float scalar_value[256], counts[256];
  double maximum, minimum;
  int k;

  maximum = minimum = histogram_array[1];
  /* counts[0] are not considered */
  for (k=1;k<256;k++){
    scalar_value[k] = (float) k-1;
    counts[k] = (float) histogram_array[k-1];
    minimum = ((double)counts[k] < minimum) ? (double)counts[k] : minimum;
    maximum = ((double)counts[k] > maximum) ? (double)counts[k] : maximum;
  }
  histogram_object = create_histogram_object();
  fl_set_xyplot_data(histogram_object->histogram_xyplot, scalar_value, counts, 255, "","","");
  fl_set_xyplot_xtics(histogram_object->histogram_xyplot, 14, 3);
  fl_set_xyplot_ytics(histogram_object->histogram_xyplot, 14, 4);
  fl_set_xyplot_ybounds(histogram_object->histogram_xyplot, 
                        (double)histogram_min, (double)histogram_max);
  /* show the form */
  fl_show_form(histogram_object->histogram_form,FL_PLACE_MOUSE,FL_TRANSIENT,"Histogram");
  fl_check_forms();
}
// ==================================================================

void selectClassifiedFile(const char *filename)
/************************************************************************/
/************************************************************************/
{
  fl_set_fselector_placement(FL_PLACE_FREE);
  filename = fl_show_fselector("Load Classified File", 0, "*.cv",0);
  if (filename != NULL){
    printf("Selected filename is now %s\n",filename);
    load_classified_volume((char *)filename);
    multi_rotate_image(1,0.0,1);
  }
  else{
    show_message("No file name selected \n","","");
  }
}

// ==================================================================

void selectUnclassifiedFile(const char *filename)
/************************************************************************/
/************************************************************************/
{
  fl_set_fselector_placement(FL_PLACE_FREE);
  filename = fl_show_fselector("Load Unclassified File", 0, "*.dat",0);
  if (filename != NULL){
    printf("Selected filename is now %s\n",filename);
    load_and_classify_volume((char *)filename);
    multi_rotate_image(1,0.0,1);
  }
  else{
    show_message("No file name selected \n","","");
  }
}

// ==================================================================
void storeImage(const char *filename)
/************************************************************************/
/************************************************************************/
{
  char filename_type [6];
  // At first only the extension has been stored in the filename (pgm or ppm)
  sprintf(filename_type,"*.%s",filename);
  fl_set_fselector_placement(FL_PLACE_FREE);
  filename = fl_show_fselector("Store Image File", 0, filename_type,0);
  if (filename != NULL){
    printf("Selected image filename is %s\n",filename);
    StorePGM((char *)filename);
  }
  else{
    show_message("No image file name selected \n","","");
  }

}
// ==================================================================
void selectStoreSequence(const char *filename)
/************************************************************************/
/************************************************************************/
{
  filename = fl_show_fselector("Select Sequence Name", 0, "*.0",0);
  if (filename != NULL){
    set_image_sequence_name((char *)filename);
  }
  else{
    show_message("No sequence name selected \n","","");
  }
}
// ==================================================================
void selectClassifiedLoadSequence(const char *filename)
/************************************************************************/
/************************************************************************/
{
  filename = fl_show_fselector("Select Sequence Name", 0, "*.cv.0",0);
  if (filename != NULL){
    load_and_render_sequence((char *)filename);
  }
  else{
    show_message("No sequence name selected \n","","");
  }
}
// ==================================================================
void selectUnclassifiedLoadSequence(const char *filename)
/************************************************************************/
/************************************************************************/
{
  filename = fl_show_fselector("Select Sequence Name", 0, "*.dat.0",0);
  if (filename != NULL){
    load_unclassified_sequence((char *)filename);
  }
  else{
    show_message("No sequence name selected \n","","");
  }
}
