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
   file name .................  render.C
   author ....................  Roland Niemeier
   version ...................  0.9beta0
   date of last change .......  25/07/97
   description : interface for direct volume visualization of scalar fields. 
                 The application development is intended for 
                 the interactive visualization of volume data sets. 
                 It is based on Open Inventor, VolPack, Xforms.
   purpose of render.C: Collection of routines that use volpack routines 

*/

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <iostream.h>
#include <new.h>
#include <bstring.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>

#include "guiforms.h"
#include "ImgRenderAction.h"
#include "global.h"
#include "viewer.h"
#include "render.h"

// volpack specific headers and data definition
#include <vp_global.h>
#include <volpack.h>

// definition of variables for local (routines in this file) and extern usage

/* rendering and drawing is only usefull if a volume is loaded */
extern Boolean DrawFlag=0;

/* parameters for rotation around an axis */
extern int number_of_rotations=9; 
extern double rotation_increment=10.0; 

/* parameters for histogram boundaries (y-axis) */
extern int histogram_minimum=0.0; 
extern int histogram_maximum=40.0; 

/* rendering context */
vpContext *vpc;	

/* either grayvalue or rgb images */ 
const int COLOR_CHANNELS = 3;
enum ColorModes {
LUMINANCE_MODE = 0,
RGB_MODE = 1
};
ColorModes CurrentColorMode=LUMINANCE_MODE;

/* definition of number of ramp points for classification ramps */
extern int n_scalar=5, n_gradient=5;

// declaration of extern variables for local usage

/* pointer to viewer, etc. as defined in global.h */
extern GlobalInfoTyp global;
extern unsigned char FORCEDRAW;


// definition of local variables

typedef struct {		/* contents of a voxel */
   short normal;		/*   encoded surface normal vector */
   unsigned char density;	/*   original density */
   unsigned char gradient;	/*   original gradient */
} RawVoxel;

RawVoxel *dummy_voxel;

#define BYTES_PER_VOXEL	sizeof(RawVoxel)	/* voxel size in bytes */
#define VOXEL_FIELDS	3	/* number of fields in voxel */
#define SHADE_FIELDS	2	/* number of fields used for shading
				   (normal and density); must be the
				   1st fields of RawVoxel */
#define CLSFY_FIELDS	2	/* number of fields used for classifying
				   (density and gradient); can be any fields
				   in the RawVoxel */

/* first voxel field (encoded normalized gradient vector) */
#define NORMAL_FIELD	0
#define NORMAL_OFFSET	vpFieldOffset(dummy_voxel, normal)
#define NORMAL_SIZE	sizeof(short)
#define NORMAL_MAX	VP_NORM_MAX

/* second voxel field of raw voxels (density) */
#define DENSITY_FIELD	1
#define DENSITY_OFFSET	vpFieldOffset(dummy_voxel, density)
#define DENSITY_SIZE	sizeof(unsigned char)
#define DENSITY_MAX	255

/* third voxel field of raw voxels (gradient magnitude) */
#define GRADIENT_FIELD	2
#define GRADIENT_OFFSET	vpFieldOffset(dummy_voxel, gradient)
#define GRADIENT_SIZE	sizeof(unsigned char)
#define GRADIENT_MAX	VP_GRAD_MAX

/* classification parameters */
#define DENSITY_PARAM		0		
#define GRADIENT_PARAM		1

/* octree parameters for fast classification */
#define OCTREE_DENSITY_THRESH	4
#define OCTREE_GRADIENT_THRESH	4
#define OCTREE_BASE_NODE_SIZE	4


/* file descriptors for raw volume, octree and classified volume */
int volume_fd, octree_fd, clvolume_fd;	

/* booleans for usage of octree data, density data or classified volume data */
int use_octree=0, use_original_data = 1, use_clvolume = 1;	

/* classification specific variables and arrays */
/* arrays for the classification ramps */
/* opacity as a function of density */
int DensityRampX[] =    {  0, 64, 128, 191, 255};
float DensityRampY[] =  {0.0, 0.25, 0.5, 0.75, 1.0};
float density_ramp[DENSITY_MAX+1];	
/* opacity as a function of gradient magnitude */
int GradientRampX[] =   {  0, 10, 20, 200, 221};
float GradientRampY[] = {1.0, 1.0, 1.0, 1.0, 1.0};
float gradient_ramp[GRADIENT_MAX+1];	
/* arrays for modifying the ramps with active xyplots */
float SRX[256], SRY[256], GRX[256], GRY[256];

/* material and light specific constants, variables or arrays */
static int current_material = 1; // for modification of material properties
static int current_light = 1; // for modification of light properties 
const int MAX_MATERIALS = 4; // Maximum number of materials that can be used 
// Currently restricted to 4, 
// see also functions modify_material_properties and 
// modify_material_numbers that must be adapted if this value is changed
int NUM_MATERIALS = 1; // number of materials that are used 
static int VP_MATERIAL; // VolPack material counts start from 0 
enum EditorModes {
INVENTOR_EDITOR = 0,
FORMS_EDITOR = 1
};
static EditorModes light_editor=FORMS_EDITOR;

/* initial_points are used for initialization of material weight ramps */
const int initial_points = 10; 
/* bookkeeper for number of ramp points for each material weight ramp */
static int material_weight_points[MAX_MATERIALS];
float WeightPointsX[MAX_MATERIALS][VP_SCALAR_MAX+1],
      WeightPointsY[MAX_MATERIALS][VP_SCALAR_MAX+1];
float WX[VP_SCALAR_MAX+1]={0.,10.,20.,30.,50.,100.,150.,200.,210.,VP_SCALAR_MAX};
float WY[VP_SCALAR_MAX+1]={0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
/* materials weight as a function of density */
float weight_ramp[VP_SCALAR_MAX+1]; 


/* Lookup tables used by VolPack's Shader
 * a test showed that the deallocation and new allocation leads to more
 * time consuming during rendering outweighting memory economics!?
*/
/* shading lookup table for grayscale mode */
float shade_table[(NORMAL_MAX+1)*MAX_MATERIALS];			
int shade_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*sizeof(float);
/* shading lookup table for rgb mode */
float color_shade_table[(NORMAL_MAX+1)*MAX_MATERIALS*COLOR_CHANNELS];   
int color_shade_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*COLOR_CHANNELS*sizeof(float);
/* materials weights lookup table */
float weight_table[(VP_SCALAR_MAX+1)*MAX_MATERIALS];
int weight_table_size = (VP_SCALAR_MAX+1)*NUM_MATERIALS*sizeof(float);
/* shadow lookup table for grayscale mode, currently not used */
//float shadow_table[(NORMAL_MAX+1)*MAX_MATERIALS];
//int shadow_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*sizeof(float);
/* shadow lookup table for rgb mode, currently not used */
// float color_shadow_table[(NORMAL_MAX+1)*MAX_MATERIALS*COLOR_CHANNELS];
// int color_shadow_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*COLOR_CHANNELS*sizeof(float);

/* image array and image size parameters */
const int IMAGE_WIDTH = 1280;	// maximum image width                  
const int IMAGE_HEIGHT = 1024;	// maximum image height
/* image buffer ([2] for stereo mode) */
static unsigned char image[2][IMAGE_WIDTH*IMAGE_HEIGHT*COLOR_CHANNELS];	
/* actual image width and height */
static int image_width, image_height;

/* volume sizes */
static int data_xlen = 128, data_ylen = 128, data_zlen = 128, data_header = 0;
/* a typical density file 128^3 data bytes, no header is assumed for 
* initialization (this is variable, the user is prompted if using a 
* new unclassified volume)
*/

/*  file name variables */
static char current_filename[256], // classified volume file name *.cv 
            density_filename[256], // density file name *.dat 
         raw_volume_filename[256], // unclassified volume file name *.rv 
             octree_filename[256], // octree file name *.oct
	  sequence_prototype[256]; // file name for sequences *.0
char output_filename[256]; // file name for image storage *.ppm, *.pgm

/* variables for the storage of images */
int sequence_number=0; // appendix of sequence_prototype
static unsigned char BoolSaveImages = 0; // if true then save the image 

/* relative sizes of a voxel */
static double scale_x=1.0, scale_y=1.0, scale_z=1.0;

/* inventor light editor view */
double iv_light_x, iv_light_y, iv_light_z;

//  declaration of functions defined in this file is done by include render.h 
//void StorePGM(char *output_filename);
//void new_classification();
//void multi_rotate_image(int axis, float angle, int number_of_rotations);
//unsigned char *rotate_image(int axis, float angle);


// definition of functions
void callback_shader(void *voxel, float *value)
/************************************************************************/
/* Callback function for shading					*/
/************************************************************************/
{
// this is just a trivial example for callback shading using
// the density value multiplied by a constant as color value
 *value = ByteField(voxel,DENSITY_OFFSET)*2.0;
}


void initialize_volpack()
/************************************************************************/
/* Initialization for volume rendering with volpack			*/
/************************************************************************/
{
    printf("Initializing volpack ...\n");
    /* create a context */
    vpc = vpCreateContext();

    /* don't use octree for classification */
    use_octree = 0;

    /* initialize the materials */
    /* for starting we only initialize two material weights specifically */
    for (int i=0;i<VP_SCALAR_MAX; i=i+2){
      weight_table[i] = 1.0;
      weight_table[i+1] = 0.0;
    }
    for (i=VP_SCALAR_MAX;i<=2*VP_SCALAR_MAX; i=i+2){
      weight_table[i] = 1.0;
      weight_table[i+1] = 1.0;
    }
    material_weight_points[0]=initial_points;
    material_weight_points[1]=initial_points;
    for(int k = 0;k<initial_points/5;k++){
      WeightPointsX[0][k] = VP_SCALAR_MAX*((float)k/(float)(initial_points-1));
      WeightPointsY[0][k] = 1.0;
      WeightPointsX[1][k] = VP_SCALAR_MAX*((float)k/(float)(initial_points-1));
      WeightPointsY[1][k] = 0.0;
    }
    for(k = initial_points/5;k<initial_points;k++){
      WeightPointsX[0][k] = VP_SCALAR_MAX*((float)k/(float)(initial_points-1));
      WeightPointsY[0][k] = 1.0;
      WeightPointsX[1][k] = VP_SCALAR_MAX*((float)k/(float)(initial_points-1));
      WeightPointsY[1][k] = 1.0;
    }
    // initialize the other materials
    // WeightPointsX[MAX_MATERIALS][VP_SCALAR_MAX+1],
    // and WeightPointsY[MAX_MATERIALS][VP_SCALAR_MAX+1];
    for(i=2; i<MAX_MATERIALS; i++){
      material_weight_points[i]=initial_points;
      for(int k = 0;k<initial_points;k++){
        WeightPointsX[i][k] = VP_SCALAR_MAX*((float)k/(float)(initial_points-1));
        WeightPointsY[i][k] = (float)k/(float)(initial_points-1);
      }
    }
    

    /* initialize the scalar volume array */
    vpSetVolumeSize(vpc, data_xlen, data_ylen, data_zlen);

    /* initialize the voxel size */
    vpSetVoxelSize(vpc, BYTES_PER_VOXEL, VOXEL_FIELDS,
		   SHADE_FIELDS, CLSFY_FIELDS);

    /* initialize the voxel fields */
    vpSetVoxelField(vpc, NORMAL_FIELD, NORMAL_SIZE, NORMAL_OFFSET,
		    NORMAL_MAX);
    vpSetVoxelField(vpc, DENSITY_FIELD, DENSITY_SIZE, DENSITY_OFFSET,
		    DENSITY_MAX);
    vpSetVoxelField(vpc, GRADIENT_FIELD, GRADIENT_SIZE, GRADIENT_OFFSET,
		    GRADIENT_MAX);

    /* set the classification function */
    vpRamp(density_ramp, sizeof(float), n_scalar, DensityRampX,
	   DensityRampY);
    vpSetClassifierTable(vpc, DENSITY_PARAM, DENSITY_FIELD, density_ramp,
			 sizeof(density_ramp));
    vpRamp(gradient_ramp, sizeof(float), n_gradient,
	   GradientRampX, GradientRampY);
    vpSetClassifierTable(vpc, GRADIENT_PARAM, GRADIENT_FIELD,
			 gradient_ramp, sizeof(gradient_ramp));

    /* voxels with an opacity lower than VP_MIN_VOXEL_OPACITY 
       are not used for compositing */
    vpSetd(vpc, VP_MIN_VOXEL_OPACITY, 0.05);

    /* set the shading lookup table */
    switch (CurrentColorMode){

      case LUMINANCE_MODE: vpSetLookupShader(vpc, 1, NUM_MATERIALS, 
                      NORMAL_FIELD, shade_table, shade_table_size, 
                      DENSITY_FIELD, weight_table, weight_table_size); 
                      break;
      case RGB_MODE: vpSetLookupShader(vpc, COLOR_CHANNELS, NUM_MATERIALS, 
                      NORMAL_FIELD, color_shade_table, color_shade_table_size, 
                      DENSITY_FIELD, weight_table, weight_table_size); 
                      break;
      default: fprintf(stderr,"Error: Unknown CurrentColorMode");
                      exit(1);
    }

    /* initialize the material parameters */
    vpSetMaterial(vpc, VP_MATERIAL0, VP_AMBIENT, VP_BOTH_SIDES,0.39, 0.0, 0.0);
    vpSetMaterial(vpc, VP_MATERIAL0, VP_DIFFUSE, VP_BOTH_SIDES,0.32, 0.0, 0.0);
    vpSetMaterial(vpc, VP_MATERIAL0, VP_SPECULAR, VP_BOTH_SIDES,0.37, 0.0, 0.0);
    vpSetMaterial(vpc, VP_MATERIAL0, VP_SHINYNESS, VP_BOTH_SIDES,4.0,0.0,0.0);
//    vpSetMaterial(vpc, VP_MATERIAL0, VP_AMBIENT, VP_BOTH_SIDES,0.0, 0.0, 0.0);
//    vpSetMaterial(vpc, VP_MATERIAL0, VP_DIFFUSE, VP_BOTH_SIDES,0.0, 0.0, 0.0);
//    vpSetMaterial(vpc, VP_MATERIAL0, VP_SPECULAR, VP_BOTH_SIDES,0.8, 0.0, 0.0);
//    vpSetMaterial(vpc, VP_MATERIAL0, VP_SHINYNESS, VP_BOTH_SIDES,4.0,0.0,0.0);

    vpSetMaterial(vpc, VP_MATERIAL3, VP_AMBIENT, VP_BOTH_SIDES,0.0, 0.45, 0.0);
    vpSetMaterial(vpc, VP_MATERIAL3, VP_DIFFUSE, VP_BOTH_SIDES,0.0, 0.45, 0.0);
    vpSetMaterial(vpc, VP_MATERIAL3, VP_SPECULAR, VP_BOTH_SIDES,0.0, 0.7, 0.0);
    vpSetMaterial(vpc, VP_MATERIAL3, VP_SHINYNESS, VP_BOTH_SIDES,4.0,0.0,0.0);

    vpSetMaterial(vpc, VP_MATERIAL2, VP_AMBIENT, VP_BOTH_SIDES,0.0, 0.0, 0.1);
    vpSetMaterial(vpc, VP_MATERIAL2, VP_DIFFUSE, VP_BOTH_SIDES,0.0, 0.0, 0.15);
    vpSetMaterial(vpc, VP_MATERIAL2, VP_SPECULAR, VP_BOTH_SIDES,0.0, 0.0, 0.6);
    vpSetMaterial(vpc, VP_MATERIAL2, VP_SHINYNESS, VP_BOTH_SIDES,4.0,0.0,0.0);

    vpSetMaterial(vpc, VP_MATERIAL1, VP_AMBIENT, VP_BOTH_SIDES,0.1, 0.1, 0.1);
    vpSetMaterial(vpc, VP_MATERIAL1, VP_DIFFUSE, VP_BOTH_SIDES,0.15, 0.15, 0.15);
    vpSetMaterial(vpc, VP_MATERIAL1, VP_SPECULAR, VP_BOTH_SIDES,0.6, 0.6, 0.6);
    vpSetMaterial(vpc, VP_MATERIAL1, VP_SHINYNESS, VP_BOTH_SIDES,4.0,0.0,0.0);

    /* initialize one light  */
    vpSetLight(vpc, VP_LIGHT0, VP_DIRECTION, 0.3, 0.3, 1.0);
    iv_light_x = 0.3;
    iv_light_y = 0.3;
    iv_light_z = 1.0;
    vpSetLight(vpc, VP_LIGHT0, VP_COLOR, 1.0, 1.0, 1.0);
    vpEnable(vpc, VP_LIGHT0, 1);

    /* set the initial viewing parameters */
    DrawFlag = 0;
    vpSeti(vpc, VP_CONCAT_MODE, VP_CONCAT_LEFT);
    vpCurrentMatrix(vpc, VP_PROJECT);
    vpIdentityMatrix(vpc);
    vpWindow(vpc, VP_PARALLEL, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5);
    vpCurrentMatrix(vpc, VP_MODEL);
    vpIdentityMatrix(vpc);
    vpCurrentMatrix(vpc, VP_VIEW);
    vpIdentityMatrix(vpc);
    vpCurrentMatrix(vpc, VP_MODEL);

    /* set the image buffer */
    switch (CurrentColorMode){
      case LUMINANCE_MODE: vpSetImage(vpc, image[0], IMAGE_WIDTH, IMAGE_HEIGHT,
 					  IMAGE_WIDTH, VP_LUMINANCE);
                      break;
      case RGB_MODE: vpSetImage(vpc, image[0], IMAGE_WIDTH, IMAGE_HEIGHT,
	      				 IMAGE_WIDTH*3, VP_RGB);
                      break;
      default: fprintf(stderr,"Error: Unknown CurrentColorMode");
                      exit(1);
    }

    /* set the maximum ray opacity */
    vpSetd(vpc, VP_MAX_RAY_OPACITY, 0.95);
}
// ---------------------------------------------------------------

void multi_rotate_image(int axis, float angle, int number_of_rotations)
/************************************************************************/
/* rotates around main axis by calling rotate_image			*/
/************************************************************************/
{
  unsigned char *image_pointer = image[0];
  if (DrawFlag){
    for(int i=0;i<number_of_rotations;i++){
      /* render with volpack, draw with GL */
      global.viewer->getImgRenderAction()->RotateImage(axis,-angle,image_pointer);
      /* store image if storage is active */
      if (BoolSaveImages){
        sprintf(output_filename, "%s.%d", sequence_prototype,sequence_number++);
        StorePGM(output_filename);
      }
    }
  }
}
// ---------------------------------------------------------------

void set_image_sequence_name(char *sequence_filename)
/************************************************************************/
/* set the name for the image sequence		 			*/
/************************************************************************/
{
  char *dot_position;
  char filename[256];
  sprintf(sequence_prototype,"%s",sequence_filename);
  strcpy(filename,sequence_prototype);
  dot_position = strstr(filename,".");
  *dot_position = '\0';
  sprintf(sequence_prototype,"%s",filename);
  sequence_number = 0;
}
// ---------------------------------------------------------------

void load_and_render_sequence(char *sequence_filename)
/************************************************************************/
/* load and render a sequence of classified volumes			*/
/************************************************************************/
{
  char *last_dot_position, *current_position;
  char filename[256], prototype[256];
  int i=0;
  FILE *fp;
  char file_exists=1;
  sprintf(prototype,"%s",sequence_filename);
  current_position = (char *)prototype;
  while (current_position != NULL){
    last_dot_position = current_position;
    // caution: last position of sequence filenname should not be a dot
    current_position++;
    current_position = strstr(current_position,".");
    if (current_position == NULL || current_position == last_dot_position){
      *last_dot_position = '\0';
    }
  }
  printf("current prototype is %s\n",prototype);
  while(file_exists){
    sprintf(filename,"%s.%d",prototype,i);
    if ((fp = fopen(filename, "r")) == NULL){
      file_exists = 0;
      i--;
    }
    else{
      fclose(fp);
      load_classified_volume(filename);
      FORCEDRAW = 1;
      multi_rotate_image(1,0.0,1);
      /* store image if storage is active */
      if (BoolSaveImages){
        printf("filename is now %s\n",filename);
        sprintf(output_filename, "%s.img", filename);
        StorePGM(output_filename);
      }
      i++;
    }
  }
  //update filenames
  sprintf(current_filename,"%s",current_filename);
  sprintf(octree_filename,"%s.%d",octree_filename,i);
  sprintf(raw_volume_filename,"%s.%d",raw_volume_filename,i);
  sprintf(density_filename,"%s.%d",density_filename,i);  
  printf("current_filename is now %s\n",current_filename);
  printf("octree_filename is now %s\n",octree_filename);
  printf("raw_volume_filename is now %s\n",raw_volume_filename);
  printf("density_filename is now %s\n",density_filename);
}
// ---------------------------------------------------------------

void load_unclassified_sequence(char *sequence_filename)
/************************************************************************/
/* load and render a sequence of classified volumes			*/
/************************************************************************/
{
  char *first_dot_position, *last_dot_position, *current_position;
  char filename[256], prototype[256];
  int i=0;
  FILE *fp;
  char file_exists=1;
  sprintf(prototype,"%s",sequence_filename);
  current_position = (char *)prototype;
  while (current_position != NULL){
    last_dot_position = current_position;
    // caution: last position of sequence filenname should not be a dot
    current_position++;
    current_position = strstr(current_position,".");
    if (current_position == NULL || current_position == last_dot_position){
      *last_dot_position = '\0';
    }
  }
  printf("current prototype is %s\n",prototype);
  while(file_exists){
    sprintf(filename,"%s.%d",prototype,i);
    if ((fp = fopen(filename, "r")) == NULL){
      file_exists = 0;
      i--;
    }
    else{
      fclose(fp);
      //update filenames
      sprintf(density_filename,"%s",filename);  
      first_dot_position = strstr(filename,".");
      *first_dot_position = '\0';
      sprintf(current_filename,"%s.cv.%d",filename,i);
      sprintf(octree_filename,"%s.oct.%d",filename,i);
      sprintf(raw_volume_filename,"%s.rv.%d",filename,i);
      printf("current_filename is now %s\n",current_filename);
      printf("octree_filename is now %s\n",octree_filename);
      printf("raw_volume_filename is now %s\n",raw_volume_filename);
      printf("density_filename is now %s\n",density_filename);
      load_and_classify_next_volume(density_filename);
      //new_classification();
      FORCEDRAW = 1;
      multi_rotate_image(1,0.0,1);
      /* store image if storage is active */
      if (BoolSaveImages){
        printf("filename is now %s\n",filename);
        sprintf(output_filename, "%s.img", filename);
        StorePGM(output_filename);
      }
      i++;
    }
  }
}
// ---------------------------------------------------------------

void build_octree()
/************************************************************************/
/* makes an octree for fast classification				*/
/************************************************************************/
{
    int volume_fd;	/* file descriptor for volume data (input) */
    int octree_fd;	/* file descriptor for octree (output) */
    /* load the volume */
    if ((volume_fd = open(raw_volume_filename, 0)) < 0) {
        show_message("Could not open file:", raw_volume_filename, 
                     "No octree is constructed");
        return;
    }
    if (vpLoadRawVolume(vpc, volume_fd) != VP_OK) {
        show_message("VolPack error:", 
                     vpGetErrorString(vpGetError(vpc)),
                     "No octree is constructed");
	return;
    }
    close(volume_fd);
    /* compute the octree */
    /* set the classification function */
    vpSetClassifierTable(vpc, DENSITY_PARAM, DENSITY_FIELD, NULL, 0);
    vpSetClassifierTable(vpc, GRADIENT_PARAM, GRADIENT_FIELD, NULL, 0);
    vpMinMaxOctreeThreshold(vpc, DENSITY_PARAM, OCTREE_DENSITY_THRESH);
    vpMinMaxOctreeThreshold(vpc, GRADIENT_PARAM, OCTREE_GRADIENT_THRESH);
    if (vpCreateMinMaxOctree(vpc, 0, OCTREE_BASE_NODE_SIZE) != VP_OK) {
        show_message("VolPack error:",
                     vpGetErrorString(vpGetError(vpc)),
                     "No octree is created (switch use octree off!)");
        return;
    }
    /* store octree in a file */
    if ((octree_fd = creat(octree_filename, 0644)) < 0) {
        show_message("Could not open new octree file:",
                      octree_filename,
                      "");
        return;
    }
    else{
      if (vpStoreMinMaxOctree(vpc, octree_fd) != VP_OK) {
  	show_message("VolPack error: %s\n",
		"vpGetErrorString(vpGetError(vpc))",
                "switch use octree off!");
	return;
      }
      close(octree_fd);
    }
}
// ---------------------------------------------------------------

void build_raw_volume()
/************************************************************************/
/* makes an raw volume (voxels containing arrays for 			*/
/* density, gradient and normal vectors)				*/
/************************************************************************/
{
    /* describe the layout of the volume */
    vpSetVolumeSize(vpc, data_xlen, data_ylen, data_zlen);
    vpSetVoxelSize(vpc, BYTES_PER_VOXEL, VOXEL_FIELDS, SHADE_FIELDS, CLSFY_FIELDS);
    vpSetVoxelField(vpc, NORMAL_FIELD, NORMAL_SIZE, NORMAL_OFFSET, NORMAL_MAX);
    vpSetVoxelField(vpc, DENSITY_FIELD, DENSITY_SIZE, DENSITY_OFFSET, DENSITY_MAX);
    vpSetVoxelField(vpc, GRADIENT_FIELD, GRADIENT_SIZE, GRADIENT_OFFSET, GRADIENT_MAX);

    /* allocate space for the raw data and the volume */
    void *volume;
    int density_fd, volume_fd;
    unsigned density_size = data_xlen * data_ylen * data_zlen;
    unsigned char *density = new unsigned char [density_size];
    unsigned volume_size = data_xlen * data_ylen * data_zlen * BYTES_PER_VOXEL;
    volume = malloc(volume_size);
    if (density == NULL || volume == NULL) {
      show_message("Can not build raw volume: ","out of memory","");
      return;
    }
    vpSetRawVoxels(vpc, volume, volume_size, BYTES_PER_VOXEL,
		   data_xlen * BYTES_PER_VOXEL,
		   data_ylen * data_xlen * BYTES_PER_VOXEL);

    /* load the raw data */
    if ((density_fd = open(density_filename, 0)) < 0) {
      show_message("Can not open scalar data file:",
                   density_filename, 
                   "No raw volume will be constructed");
      return;
    }
    if (lseek(density_fd, data_header, 0) < 0) {
      show_message("could not skip header of file:", 
                   density_filename, 
                   "No raw volume file will be constructed!");
      return;
    }
    if (read(density_fd, density, density_size) != density_size) {
      show_message("could not read data from file:", 
                   density_filename, 
                   "No raw volume file will be constructed!");
      return;
    }
    close(density_fd);

    /* compute surface normals (for shading) and
       gradient magnitudes (for classification) */
    if (vpVolumeNormals(vpc, (unsigned char *)density, density_size, DENSITY_FIELD,
			GRADIENT_FIELD, NORMAL_FIELD) != VP_OK) {
      show_message("VolPack error:",
                   vpGetErrorString(vpGetError(vpc)),
                   "No raw volume created!");
      return;
    }

    /* store volume in a file */
    if ((volume_fd = creat(raw_volume_filename, 0644)) < 0) {
      show_message("could not create new file:", raw_volume_filename, 
      "No raw volume file is created or stored");
      return;
    }
    if (vpStoreRawVolume(vpc, volume_fd) != VP_OK) {
      show_message("VolPack error:",
                   vpGetErrorString(vpGetError(vpc)),
                   "No raw volume file is stored");
      return;
    }
    close(volume_fd);
    delete density;
}

// ---------------------------------------------------------------

void load_classified_volume(char *cl_volume_filename)
/************************************************************************/
/* load a new classified volume						*/
/************************************************************************/
{
  char *dot_position;
  char filename[256];
  int clvolume_fd; 
  if ((clvolume_fd = open(cl_volume_filename, 0)) < 0) {
    show_message("could not open classified volume from file:",
                  cl_volume_filename,"");
    DrawFlag = 0;
  }
  else{
    vpDestroyClassifiedVolume(vpc); 
    if (vpLoadClassifiedVolume(vpc, clvolume_fd) != VP_OK) {
      fprintf(stderr,"Error: Can not load the classified volume!\n");
    }
    data_xlen = vpc->xlen;
    data_ylen = vpc->ylen;
    data_zlen = vpc->zlen;
    printf("New dimensions are %d %d %d\n",data_xlen,data_ylen,data_zlen);
    sprintf(current_filename,"%s",cl_volume_filename);
    strcpy(filename,current_filename);
    dot_position = strstr(filename,".");
    *dot_position = '\0';
    sprintf(octree_filename,"%s.oct",filename);
    sprintf(raw_volume_filename,"%s.rv",filename);
    sprintf(density_filename,"%s.dat",filename);
    printf("new classified volume filename is: %s\n", current_filename);
    printf("new octree_filename is: %s\n", octree_filename);
    printf("new raw_volume_filename is: %s\n", raw_volume_filename);
    printf("new density_filename is: %s\n", density_filename);
    DrawFlag = 1;
  }
}
// ---------------------------------------------------------------

void change_filenames(char *filename)
/************************************************/
/* updates the filenames			*/
/************************************************/
{
  char *dot_position;
  int density_fd; 
  if ((density_fd = open(filename, 0, O_RDONLY)) < 0) {
    perror("open");
    fprintf(stderr, "could not open %s\n", filename);
    DrawFlag = 0;
  }
  else{ 
    strcpy(density_filename,filename);
    dot_position = strstr(filename,".");
    *dot_position = '\0';
    sprintf(octree_filename,"%s.oct",filename);
    sprintf(raw_volume_filename,"%s.rv",filename);
    sprintf(current_filename,"%s.cv",filename);
    printf("new octree_filename is: %s\n", octree_filename);
    printf("new raw_volume_filename is: %s\n", raw_volume_filename);
    printf("new current_filename is: %s\n", current_filename);
    printf("new density_filename is: %s\n", density_filename);
  }
}

// ---------------------------------------------------------------

void automatic_load_image(char *filename, int x_dim, int y_dim)
/************************************************************************/
/* load a new data set from a 2D Simulation and prepare it for display	*/
/************************************************************************/
{
  struct stat statbuf; long selected_size;
  int image_fd, i, k;
  if ((image_fd = open(filename, 0, O_RDONLY)) < 0) {
    perror("open");
    fprintf(stderr, "could not open %s\n", filename);
    DrawFlag = 0;
  }
  else{ 
    //set_dimensions(&data_header, &data_xlen, &data_ylen, &data_zlen, 
    //  &apply_button_pushed);
    if(fstat(image_fd,&statbuf)==-1){ 
      fprintf(stderr,"In subroutine automatic_load_image: fstat returns -1\n");
      exit(1);
    };
    if( (statbuf.st_mode & S_IFMT) == S_IFDIR){
      fprintf(stderr,
              "In subroutine automatic_load_image: wrong file mode of file %s\n",
              filename);
      exit(1);
    } 
    selected_size = x_dim * y_dim;
    if (selected_size <= statbuf.st_size){
      if (selected_size < statbuf.st_size){
        char message_buf[16];
        sprintf(message_buf,"%d",selected_size);
        show_message("Warning: Selected filesize:",message_buf,"is smaller then real size!\n");
      }
      DrawFlag = 1;
    }
    else{
      show_message("Error: ","Selected size larger than file size!",
          "No new classification or rendering");
      DrawFlag = 0;
    }
  }
  if (DrawFlag){
    unsigned char *image_array = new unsigned char[data_xlen*data_ylen];
    if (read(image_fd, image_array, data_xlen*data_ylen) != data_xlen*data_ylen) {
      show_message("could not read data from file:",filename,"");
      DrawFlag=0; return;
    }
    /* without scaling, centering just put it to an image edge */
     global.viewer->getImgRenderAction()->getImgSize(&image_width,&image_height);
     printf("image_width and image_height is %d %d\n",image_width,image_height);
     printf("data_xlen and data_ylen is %d %d\n",data_xlen,data_ylen);
     if (data_ylen < image_height){
      for(i=0;i<data_ylen;i++){
        if (data_xlen < image_width){
          for(k=0;k<data_xlen;k++){
            image[0][k+image_width*i]=image_array[k+data_xlen*i];
          }
          for(k=data_xlen;k<image_width;k++){
            image[0][k+image_width*i]=0;
          }
        }
        else{
          for(k=0;k<image_width;k++){
            image[0][k+image_width*i]=image_array[k+data_xlen*i];
          }
        }
      }
      for(i=data_ylen;i<image_height;i++){
        for(k=0;k<image_width;k++){
          image[0][k+image_width*i]=0;
        }
      }
    }
    else{
      for(i=0;i<image_height;i++){
        if (data_xlen < image_width){
          for(k=0;k<data_xlen;k++){
            image[0][k+image_width*i]=image_array[k+data_xlen*i];
          }
          for(k=data_xlen;k<image_width;k++){
            image[0][k+image_width*i]=0;
          }
        }
        else{
          for(k=0;k<image_width;k++){
            image[0][k+image_width*i]=image_array[k+data_xlen*i];
          }
        }
      }
    }
    free(image_array);
  }
}
// ---------------------------------------------------------------

void load_image(char *filename)
/************************************************************************/
/* load a new data set from a 2D Simulation and prepare it for display	*/
/************************************************************************/
{
  struct stat statbuf; long selected_size;
  unsigned char apply_button_pushed = 0;
  int image_fd, i, k;
  if ((image_fd = open(filename, 0, O_RDONLY)) < 0) {
    perror("open");
    fprintf(stderr, "could not open %s\n", filename);
    DrawFlag = 0;
  }
  else{ 
    set_dimensions(&data_header, &data_xlen, &data_ylen, &data_zlen, 
      &apply_button_pushed);
    if(fstat(image_fd,&statbuf)==-1){ 
      fprintf(stderr,"In subroutine automatic_load_image: fstat returns -1\n");
      exit(1);
    };
    if( (statbuf.st_mode & S_IFMT) == S_IFDIR){
      fprintf(stderr,
              "In subroutine load_image: wrong file mode of file %s\n",
              filename);
      exit(1);
    } 
      selected_size = data_xlen * data_ylen + data_header;
      if (selected_size <= statbuf.st_size && apply_button_pushed){
        if (selected_size < statbuf.st_size){
          char message_buf[16];
          sprintf(message_buf,"%d",selected_size);
          show_message("Warning: Selected filesize:",message_buf,"is smaller then real size!\n");
        }
        DrawFlag = 1;
      }
      else{
        show_message("Error: ","Selected size larger than file size!",
          "No new classification or rendering");
        DrawFlag = 0;
      }
  }
  if (DrawFlag){
    unsigned char *image_array = new unsigned char[data_xlen*data_ylen];
    if (read(image_fd, image_array, data_xlen*data_ylen) != data_xlen*data_ylen) {
      perror("read");
      fprintf(stderr, "could not read data from %s\n", filename);
      exit(1);
    }
    /* without scaling, centering just put it to an image edge */
     global.viewer->getImgRenderAction()->getImgSize(&image_width,&image_height);
     printf("image_width and image_height is %d %d\n",image_width,image_height);
     printf("data_xlen and data_ylen is %d %d\n",data_xlen,data_ylen);
     if (data_ylen < image_height){
      for(i=0;i<data_ylen;i++){
        if (data_xlen < image_width){
          for(k=0;k<data_xlen;k++){
            image[0][k+image_width*i]=image_array[k+data_xlen*i];
          }
          for(k=data_xlen;k<image_width;k++){
            image[0][k+image_width*i]=0;
          }
        }
        else{
          for(k=0;k<image_width;k++){
            image[0][k+image_width*i]=image_array[k+data_xlen*i];
          }
        }
      }
      for(i=data_ylen;i<image_height;i++){
        for(k=0;k<image_width;k++){
          image[0][k+image_width*i]=0;
        }
      }
    }
    else{
      for(i=0;i<image_height;i++){
        if (data_xlen < image_width){
          for(k=0;k<data_xlen;k++){
            image[0][k+image_width*i]=image_array[k+data_xlen*i];
          }
          for(k=data_xlen;k<image_width;k++){
            image[0][k+image_width*i]=0;
          }
        }
        else{
          for(k=0;k<image_width;k++){
            image[0][k+image_width*i]=image_array[k+data_xlen*i];
          }
        }
      }
    }
    free(image_array);
  }
}
// ---------------------------------------------------------------

unsigned char *get_image_pointer(void){
/************************************************************************/
/* returns the adress of the image					*/
/************************************************************************/
return image[0];
};

// ---------------------------------------------------------------


void automatic_load_and_classify_volume(char *filename, int xdim, int ydim, int zdim)
/*************************************************************************/
/* load a new data set and classify it without asking for the resolutions*/
/*************************************************************************/
{
  char *dot_position;
  struct stat statbuf; long selected_size;
  int volume_fd; 
  if ((volume_fd = open(filename, 0, O_RDONLY)) < 0) {
    perror("open");
    fprintf(stderr, "could not open %s\n", filename);
    DrawFlag = 0;
  }
  else{ 
    strcpy(density_filename,filename);
    dot_position = strstr(filename,".");
    *dot_position = '\0';
    sprintf(octree_filename,"%s.oct",filename);
    sprintf(raw_volume_filename,"%s.rv",filename);
    sprintf(current_filename,"%s.cv",filename);
    printf("new octree_filename is: %s\n", octree_filename);
    printf("new raw_volume_filename is: %s\n", raw_volume_filename);
    printf("new current_filename is: %s\n", current_filename);
    printf("new density_filename is: %s\n", density_filename);
    data_xlen = xdim; data_ylen = ydim; data_zlen = zdim;
    data_header = 0;
    if(fstat(volume_fd,&statbuf)==-1){ 
      fprintf(stderr,"In subroutine automatic_load_and_classify_volume: fstat returns -1\n");
      exit(1);
    };
    if( (statbuf.st_mode & S_IFMT) == S_IFDIR){
      fprintf(stderr,
              "In subroutine automatic_load_and_classify_volume: wrong file mode of file %s\n",
              filename);
      exit(1);
    } 
    selected_size = data_xlen * data_ylen * data_zlen + data_header;
    if (selected_size <= statbuf.st_size){
      if (selected_size < statbuf.st_size){
        char message_buf[16];
        sprintf(message_buf,"%d",selected_size);
        show_message("Warning: Selected filesize:",message_buf,"is smaller then real size!\n");
      }
      vpSetVolumeSize(vpc, data_xlen, data_ylen, data_zlen);
      new_classification();
      DrawFlag = 1;
    }
    else{
      show_message("Error: ","Selected size larger than file size!",
        "No new classification or rendering");
      DrawFlag = 0;
    }
  }
}


// ---------------------------------------------------------------
void load_and_classify_volume(char *filename)
/************************************************************************/
/* load a new data set and classify it					*/
/************************************************************************/
{
  char *dot_position;
  struct stat statbuf; long selected_size;
  unsigned char apply_button_pushed = 0;
  int volume_fd; 
  if ((volume_fd = open(filename, 0, O_RDONLY)) < 0) {
    perror("open");
    fprintf(stderr, "could not open %s\n", filename);
    DrawFlag = 0;
  }
  else{ 
    strcpy(density_filename,filename);
    dot_position = strstr(filename,".");
    *dot_position = '\0';
    sprintf(octree_filename,"%s.oct",filename);
    sprintf(raw_volume_filename,"%s.rv",filename);
    sprintf(current_filename,"%s.cv",filename);
    printf("new octree_filename is: %s\n", octree_filename);
    printf("new raw_volume_filename is: %s\n", raw_volume_filename);
    printf("new current_filename is: %s\n", current_filename);
    printf("new density_filename is: %s\n", density_filename);
    set_dimensions(&data_header, &data_xlen, &data_ylen, &data_zlen, 
      &apply_button_pushed);
    if(apply_button_pushed){
      if(fstat(volume_fd,&statbuf)==-1) exit(1);
      if( (statbuf.st_mode & S_IFMT) == S_IFDIR) exit(1);
      selected_size = data_xlen * data_ylen * data_zlen + data_header;
      if (selected_size <= statbuf.st_size && apply_button_pushed){
        if (selected_size < statbuf.st_size){
          fprintf(stderr,
          "Warning: Selected filesize %d is smaller then real size!\n",
               selected_size);
        }
        vpSetVolumeSize(vpc, data_xlen, data_ylen, data_zlen);
        new_classification();
        DrawFlag = 1;
      }
      else{
        show_message("Error:","Selected size larger than file size!","");
        DrawFlag = 0;
      }
    }
  }
}
// ---------------------------------------------------------------
void load_and_classify_next_volume(char *filename)
/************************************************************************/
/* load next data set and classify it					*/
/************************************************************************/
{
  int volume_fd; 
  if ((volume_fd = open(filename, 0, O_RDONLY)) < 0) {
    show_message("could not open classified volume file:",filename,"");
    DrawFlag = 0;
  }
  else{ 
    close(volume_fd); /* FG - we cannot leave all those files open... */
    vpSetVolumeSize(vpc, data_xlen, data_ylen, data_zlen);
    new_classification();
    DrawFlag = 1;
  }
}
// ---------------------------------------------------------------

void switch_octree_mode(void)
/************************************************************************/
/* toggles use of octree for the classification				*/
/************************************************************************/
{ 
  int fd; 
  char string1[300];
  use_octree = (use_octree) ? 0 : 1;
  use_original_data = (use_octree) ? 0 : 1;
  if (use_octree){
    // check if files are selected
    if ((fd = open(density_filename, 0)) < 0) {
      show_message("No files are currently selected!!","","");
    }
    else{
      close(fd);
      // check if the octree for the selected file exists
      if ((fd = open(octree_filename, 0)) < 0) {
      sprintf(string1,"No necessary octree file:\n %s",current_filename);
      sprintf(string1,"%s\nfor fast classification exist",string1);
      sprintf(string1,"%s\nBuild it now?",string1);
      if (show_question(string1)){
        build_raw_volume(); build_octree();
      }
    }
  else{
    close(fd);
  }
    }
  }
}
// ---------------------------------------------------------------

void switch_color_mode(void)
/************************************************************************/
/* toggles RGB - Luminance Mode			 			*/
/************************************************************************/
{
  printf("entering switch_color_mode\n");
     if (CurrentColorMode==LUMINANCE_MODE) {
       CurrentColorMode = RGB_MODE;
       vpSetLookupShader(vpc, COLOR_CHANNELS, NUM_MATERIALS, 
                      NORMAL_FIELD, color_shade_table, color_shade_table_size, 
                      DENSITY_FIELD, weight_table, weight_table_size); 
       vpSetImage(vpc, image[0], image_width, image_height, 
                       image_width*COLOR_CHANNELS, VP_RGB);
     } 
     else {
       CurrentColorMode = LUMINANCE_MODE;
       vpSetLookupShader(vpc, 1, NUM_MATERIALS, 
                      NORMAL_FIELD, shade_table, shade_table_size, 
                      DENSITY_FIELD, weight_table, weight_table_size); 
       vpSetImage(vpc, image[0], image_width, image_height, 
                       image_width, VP_LUMINANCE); 
     }
  multi_rotate_image(1,0.0,1);
}
// ---------------------------------------------------------------

void switch_image_store_mode(void)
/************************************************************************/
/* toggles storage mode of displayed images (during rotation only) 	*/
/************************************************************************/
{
  BoolSaveImages = (BoolSaveImages) ? 0 : 1;
}
// ---------------------------------------------------------------

void switch_light_editor(void)
/************************************************************************/
/* toggles storage mode of displayed images (during rotation only) 	*/
/************************************************************************/
{
  light_editor = (light_editor == FORMS_EDITOR) ? INVENTOR_EDITOR : FORMS_EDITOR;
}
// ---------------------------------------------------------------

void vpToggle(int vp_option)
/************************************************************************/
/* toggles vpEnable options						*/
/************************************************************************/
{
  int iptr, value;
  if (vpGeti(vpc, vp_option, &iptr) != VP_OK) {
    show_message("VolPack error:",vpGetErrorString(vpGetError(vpc)),"");
  }
  else{
    value = 0;
    if (iptr == 0) value = 1;
    vpEnable(vpc, vp_option, value);
  };
}
// ---------------------------------------------------------------

void set_material_ramp(int number_of_points){
/************************************************************************/
/* set the new weight_table if apply button is pushed in 		*/
/* submenu material properties						*/
/************************************************************************/
  int i,k;
  int WeightRampX[VP_SCALAR_MAX+1];
  for (i=0;i<number_of_points;i++){
    WeightRampX[i] = (int) WX[i];
  }
  material_weight_points[current_material-1] = number_of_points;
  vpRamp(weight_ramp, sizeof(float), number_of_points,
	   WeightRampX, WY);
  //for (i=0;i<=VP_SCALAR_MAX;i++){
  //  printf("weight_ramp[%d] is %f \n",i,weight_ramp[i]);
  //}
  for (k=0,i=(current_material-1);k<=VP_SCALAR_MAX;k++,i=i+NUM_MATERIALS){
    weight_table[i] = weight_ramp[k];
  }
}
// ---------------------------------------------------------------

void set_material_ramps(){
/************************************************************************/
/* set the new weight_tables if number of materials has changed 	*/
/* can be done using submenu material numbers				*/
/************************************************************************/
  int i,k,l;
  int WeightRampX[VP_SCALAR_MAX+1];
  for (l=1;l<=NUM_MATERIALS;l++){
    for (i=0;i<material_weight_points[l-1];i++){
      WeightRampX[i] = (int) WeightPointsX[l-1][i];
    }
    vpRamp(weight_ramp, sizeof(float), material_weight_points[l-1],
  	   WeightRampX, WeightPointsY[l-1]);
    for (k=0,i=(l-1);k<=VP_SCALAR_MAX;k++,i=i+NUM_MATERIALS){
      weight_table[i] = weight_ramp[k];
    }
  }
}
// ---------------------------------------------------------------

void modify_material_properties()
/************************************************************************/
/* popups the menu for modifying the material properties		*/
/************************************************************************/
{
double r_ambient, g_ambient, b_ambient, r_diffuse, g_diffuse, b_diffuse,
       r_specular, g_specular, b_specular, shinyness, dummy_shinyness;
int i;
int number_of_materials = NUM_MATERIALS;
switch (current_material){
  case 1 : VP_MATERIAL = VP_MATERIAL0;
           break;
  case 2 : VP_MATERIAL = VP_MATERIAL1;
           break;
  case 3 : VP_MATERIAL = VP_MATERIAL2;
           break;
  case 4 : VP_MATERIAL = VP_MATERIAL3;
           break;
  case 5 : VP_MATERIAL = VP_MATERIAL4;
           show_message("Sorry in this implementation",
                        "the maximum number of materials is 4","");
           return;
           //break;
  case 6 : VP_MATERIAL = VP_MATERIAL5;
           show_message("Sorry in this implementation",
                        "the maximum number of materials is 4","");
           return;
           //break;
  default: fprintf(stderr, "material number %d is not valid",current_material);
           exit(1);  
};
vpGetMaterial(vpc, VP_MATERIAL, VP_AMBIENT, VP_EXTERIOR, 
                   &r_ambient, &g_ambient, &b_ambient);
vpGetMaterial(vpc, VP_MATERIAL, VP_DIFFUSE, VP_EXTERIOR, 
                   &r_diffuse, &g_diffuse, &b_diffuse);
vpGetMaterial(vpc, VP_MATERIAL, VP_SPECULAR, VP_EXTERIOR, 
                   &r_specular, &g_specular, &b_specular);
vpGetMaterial(vpc, VP_MATERIAL, VP_SHINYNESS, VP_EXTERIOR, 
                   &shinyness, &dummy_shinyness, &dummy_shinyness);
// build WX, WY
for (i=0;i<material_weight_points[current_material-1];i++){
  WX[i] = WeightPointsX[current_material-1][i];
  WY[i] = WeightPointsY[current_material-1][i];
}
modify_material_popup(&r_ambient,&g_ambient,&b_ambient,
		      &r_diffuse,&g_diffuse,&b_diffuse,
                      &r_specular,&g_specular,&b_specular,
		      &shinyness,&current_material,&number_of_materials,
                      material_weight_points[current_material-1]);
NUM_MATERIALS = number_of_materials;

}
// ---------------------------------------------------------------

void update_material_properties(
  double *r_ambient, double *g_ambient, double *b_ambient,
  double *r_diffuse, double *g_diffuse, double *b_diffuse,
  double *r_specular, double *g_specular, double *b_specular,
  double *shinyness, int new_current_material, int *material_weights_pointer)
/************************************************************************/
/* updatess the menu for modifying the material properties		*/
/************************************************************************/
{
double dummy_shinyness;
int i;
current_material = new_current_material;
switch (current_material){
  case 1 : VP_MATERIAL = VP_MATERIAL0;
           break;
  case 2 : VP_MATERIAL = VP_MATERIAL1;
           break;
  case 3 : VP_MATERIAL = VP_MATERIAL2;
           break;
  case 4 : VP_MATERIAL = VP_MATERIAL3;
           break;
  case 5 : VP_MATERIAL = VP_MATERIAL4;
           fprintf(stderr,"Sorry in this implementation maximum of materials is 4!\n");
           return;
           //break;
  case 6 : VP_MATERIAL = VP_MATERIAL5;
           fprintf(stderr,"Sorry in this implementation maximum of materials is 4!\n");
           return;
           //break;
  default: fprintf(stderr, "material number %d is not valid",current_material);
           exit(1);  
};
vpGetMaterial(vpc, VP_MATERIAL, VP_AMBIENT, VP_EXTERIOR, 
                   r_ambient, g_ambient, b_ambient);
vpGetMaterial(vpc, VP_MATERIAL, VP_DIFFUSE, VP_EXTERIOR, 
                   r_diffuse, g_diffuse, b_diffuse);
vpGetMaterial(vpc, VP_MATERIAL, VP_SPECULAR, VP_EXTERIOR, 
                   r_specular, g_specular, b_specular);
vpGetMaterial(vpc, VP_MATERIAL, VP_SHINYNESS, VP_EXTERIOR, 
                   shinyness, &dummy_shinyness, &dummy_shinyness);
// build WX, WY
for (i=0;i<material_weight_points[current_material-1];i++){
  WX[i] = WeightPointsX[current_material-1][i];
  WY[i] = WeightPointsY[current_material-1][i];
}
*material_weights_pointer = material_weight_points[current_material-1];
printf("material_weight_points %d\n",material_weight_points[current_material-1]);
printf("number of material points is %d\n",*material_weights_pointer);
printf("current_material is nr %d\n",current_material);
}
// ---------------------------------------------------------------
void set_material_weight_points(int new_material_weight_points, int current_material){
/************************************************************************/
/************************************************************************/
  material_weight_points[current_material-1] = new_material_weight_points;
}
// ---------------------------------------------------------------

void modify_material_properties_applyCB(double r_ambient, double g_ambient, double b_ambient,
		       double r_diffuse,double g_diffuse,double b_diffuse,
                       double r_specular,double  g_specular,double b_specular,
		       double shinyness)
/************************************************************************/
/* called if apply button in modify_material_properties popup is pushed	*/
/************************************************************************/
{
  for (int i=0;i<material_weight_points[current_material-1];i++){
    WeightPointsX[current_material-1][i] = WX[i];
    WeightPointsY[current_material-1][i] = WY[i];
  }
  switch (current_material){
    case 1 : VP_MATERIAL = VP_MATERIAL0;
             break;
    case 2 : VP_MATERIAL = VP_MATERIAL1;
             break;
    case 3 : VP_MATERIAL = VP_MATERIAL2;
             break;
    case 4 : VP_MATERIAL = VP_MATERIAL3;
             break;
    case 5 : VP_MATERIAL = VP_MATERIAL4;
             fprintf(stderr,"Sorry in this implementation maximum of materials is 4!\n");
             return;
             //break;
    case 6 : VP_MATERIAL = VP_MATERIAL5;
             fprintf(stderr,"Sorry in this implementation maximum of materials is 4!\n");
             return;
             //break;
    default: fprintf(stderr, "material number %d is not valid",current_material);
             exit(1);  
  };
  vpSetMaterial(vpc, VP_MATERIAL, VP_AMBIENT, VP_BOTH_SIDES, 
                   r_ambient, g_ambient, b_ambient);
  vpSetMaterial(vpc, VP_MATERIAL, VP_DIFFUSE, VP_BOTH_SIDES, 
                   r_diffuse, g_diffuse, b_diffuse);
  vpSetMaterial(vpc, VP_MATERIAL, VP_SPECULAR, VP_BOTH_SIDES, 
                   r_specular, g_specular, b_specular);
  vpSetMaterial(vpc, VP_MATERIAL, VP_SHINYNESS, VP_BOTH_SIDES, 
                   shinyness, shinyness, shinyness);
  multi_rotate_image(1,0.0,1);
}
// ---------------------------------------------------------------

void update_material_numbers(int number_of_materials)
/********************************************************/
/* updates the total number of materials if modified    */
/* in material_properties_popupnd allocates new		*/
/********************************************************/
{
  // in the current implementation the maximum of materials is 4
  // (see the variable definition MAX_MATERIAL)
  if (number_of_materials > 4) {
    fprintf(stderr,"Sorry in this implementation maximum of materials is 4!");
    return;
  }
  // update the number of materials used
  NUM_MATERIALS = number_of_materials;
  // set the new dimension for the Lookup Shader
  shade_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*sizeof(float);
  color_shade_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*COLOR_CHANNELS*sizeof(float);
  weight_table_size = (VP_SCALAR_MAX+1)*NUM_MATERIALS*sizeof(float);
  /* first set up the weight tables */
  set_material_ramps();
  /* set the shading lookup table */
  switch (CurrentColorMode){
    case LUMINANCE_MODE: vpSetLookupShader(vpc, 1, NUM_MATERIALS, 
                      NORMAL_FIELD, shade_table, shade_table_size, 
                      DENSITY_FIELD, weight_table, weight_table_size); 
                      break;
    case RGB_MODE: vpSetLookupShader(vpc, COLOR_CHANNELS, NUM_MATERIALS, 
                      NORMAL_FIELD, color_shade_table, color_shade_table_size, 
                      DENSITY_FIELD, weight_table, weight_table_size); 
                      break;
    default: fprintf(stderr,"Error: Unknown CurrentColorMode");
                      exit(1);
  }
}
// ---------------------------------------------------------------
int get_number_of_materials(){
/****************************************************************/
/* returns the current number of materials			*/
/****************************************************************/
  return NUM_MATERIALS;
}
// ---------------------------------------------------------------

void modify_material_numbers()
/********************************************************************************/
/* calls the popup-menu for modifying the material numbers and allocates new		*/
/* memory for the lookup tables if the total number of material has changed     */
/********************************************************************************/
{
  int number_of_materials = NUM_MATERIALS;
  // calling the popup menu 
  material_numbers_sliders(&number_of_materials, &current_material);
  // in the current implementation the maximum of materials is 4
  // (see the variable definition MAX_MATERIAL)
  if (number_of_materials > 4) {
    fprintf(stderr,"Sorry in this implementation maximum of materials is 4!");
    return;
  }
  if (number_of_materials != NUM_MATERIALS){
  // update the number of materials used
    NUM_MATERIALS = number_of_materials;
  // set the new dimension for the Lookup Shader
    shade_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*sizeof(float);
    color_shade_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*COLOR_CHANNELS*sizeof(float);
    weight_table_size = (VP_SCALAR_MAX+1)*NUM_MATERIALS*sizeof(float);
  // shadowing is currently not used
  //    shadow_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*sizeof(float);
  //    color_shadow_table_size = (NORMAL_MAX+1)*NUM_MATERIALS*COLOR_CHANNELS*sizeof(float);

  // a test showed that the deallocation and new allocation leads to more
  // time consuming during rendering outweighting memory economics!?
  //    delete shade_table, color_shade_table, weight_table, shadow_table,
  //           color_shadow_table;
  //    float *shade_table = new float [(NORMAL_MAX+1)*number_of_materials];	/* shading lookup table */
  //    float *color_shade_table = new float [(NORMAL_MAX+1)*number_of_materials*COLOR_CHANNELS];	/* shading lookup table */
  //    float *weight_table = new float [(VP_SCALAR_MAX+1)*number_of_materials]; /* material weight lookup table */
  //    float *shadow_table = new float [(NORMAL_MAX+1)*number_of_materials];	/* shadow lookup table */
  //    float *color_shadow_table = new float [(NORMAL_MAX+1)*NUM_MATERIALS*COLOR_CHANNELS];	/* shadow lookup table */
  //    if ((!shade_table) || (!color_shade_table) || (!weight_table) 
  //                     || (!shadow_table) || (!color_shadow_table)){ 
  //      fprintf(stderr, "could not allocate memory for lookup tables!\n"); 
  //      exit(1);
  //    }
    /* first set up the weight tables */
    set_material_ramps();
    /* set the shading lookup table */
    switch (CurrentColorMode){
      case LUMINANCE_MODE: vpSetLookupShader(vpc, 1, NUM_MATERIALS, 
                      NORMAL_FIELD, shade_table, shade_table_size, 
                      DENSITY_FIELD, weight_table, weight_table_size); 
                      break;
      case RGB_MODE: vpSetLookupShader(vpc, COLOR_CHANNELS, NUM_MATERIALS, 
                      NORMAL_FIELD, color_shade_table, color_shade_table_size, 
                      DENSITY_FIELD, weight_table, weight_table_size); 
                      break;
      default: fprintf(stderr,"Error: Unknown CurrentColorMode");
                      exit(1);
    }
  }
}
// ---------------------------------------------------------------

void modify_light_properties()
/************************************************************************/
/* popups the menu for modifying the light properties			*/
/************************************************************************/
{
  double r_light, g_light, b_light;
  double x_direction, y_direction, z_direction;
  double x_1, y_1, z_1;
   //double view_x, view_y, view_z;
  int VP_LIGHT;
  vpMatrix4 viewing_matrix;
  switch (current_light){
    case 1 : VP_LIGHT = VP_LIGHT0;
         break;
    case 2 : VP_LIGHT = VP_LIGHT1;
         break;
    case 3 : VP_LIGHT = VP_LIGHT2;
         break;
    case 4 : VP_LIGHT = VP_LIGHT3;
         break;
    case 5 : VP_LIGHT = VP_LIGHT4;
         break;
    case 6 : VP_LIGHT = VP_LIGHT5;
         break;
    default: fprintf(stderr, "light number %d is not valid",current_light);
         exit(1);  
  };
  vpGetLight(vpc, VP_LIGHT, VP_COLOR, &r_light, &g_light, &b_light);
  vpGetLight(vpc, VP_LIGHT, VP_DIRECTION, &x_direction, 
              &y_direction, &z_direction);
  switch (light_editor){
    case INVENTOR_EDITOR:
      vpGetMatrix(vpc, VP_MODEL, viewing_matrix);
      x_1 = viewing_matrix[0][0] * x_direction +
               viewing_matrix[0][1] * y_direction +
               viewing_matrix[0][2] * z_direction;                
      y_1 = viewing_matrix[1][0] * x_direction +
               viewing_matrix[1][1] * y_direction +
               viewing_matrix[1][2] * z_direction;                
      z_1 = viewing_matrix[2][0] * x_direction +
               viewing_matrix[2][1] * y_direction +
               viewing_matrix[2][2] * z_direction;                
      iv_light_x = viewing_matrix[0][0] * x_1 +
               viewing_matrix[0][1] * y_1 +
               viewing_matrix[0][2] * z_1;                
      iv_light_y = viewing_matrix[1][0] * x_1 +
               viewing_matrix[1][1] * y_1 +
               viewing_matrix[1][2] * z_1;                
      iv_light_z = viewing_matrix[2][0] * x_1 +
               viewing_matrix[2][1] * y_1 +
               viewing_matrix[2][2] * z_1;                
      ModifyLightCB(r_light, g_light, b_light, 
              iv_light_x, iv_light_y, iv_light_z);
      break;
    case FORMS_EDITOR:
      modify_lighting_popup(&r_light, &g_light, &b_light, 
              &x_direction, &y_direction, &z_direction, current_light);
      //vpSetLight(vpc, VP_LIGHT, VP_COLOR, r_light, g_light, b_light);
      //vpSetLight(vpc, VP_LIGHT, VP_DIRECTION, x_direction, 
      //      y_direction, z_direction);
      //vpShadeTable(vpc);
      break;
  }  
}
// ---------------------------------------------------------------

void setLightColors(double r_light, double g_light, double b_light)
/************************************************************************/
/* set the colors of the current light 				*/
/************************************************************************/
{
  int VP_LIGHT;
  switch (current_light){
    case 1 : VP_LIGHT = VP_LIGHT0;
         break;
    case 2 : VP_LIGHT = VP_LIGHT1;
         break;
    case 3 : VP_LIGHT = VP_LIGHT2;
         break;
    case 4 : VP_LIGHT = VP_LIGHT3;
         break;
    case 5 : VP_LIGHT = VP_LIGHT4;
         break;
    case 6 : VP_LIGHT = VP_LIGHT5;
         break;
    default: fprintf(stderr, "light number %d is not valid",current_light);
         exit(1);  
  };
  vpSetLight(vpc, VP_LIGHT, VP_COLOR, r_light, g_light, b_light);
  vpShadeTable(vpc);
}

// ---------------------------------------------------------------

void setLightDirection(double x_direction,double  y_direction, double z_direction)
/************************************************************************/
/* set the direction of the current light 				*/
/************************************************************************/
{
  int VP_LIGHT;
  switch (current_light){
    case 1 : VP_LIGHT = VP_LIGHT0;
         break;
    case 2 : VP_LIGHT = VP_LIGHT1;
         break;
    case 3 : VP_LIGHT = VP_LIGHT2;
         break;
    case 4 : VP_LIGHT = VP_LIGHT3;
         break;
    case 5 : VP_LIGHT = VP_LIGHT4;
         break;
    case 6 : VP_LIGHT = VP_LIGHT5;
         break;
    default: fprintf(stderr, "light number %d is not valid",current_light);
         exit(1);  
  };
  vpSetLight(vpc, VP_LIGHT, VP_DIRECTION, x_direction, y_direction, z_direction);
  vpShadeTable(vpc);
}

// ---------------------------------------------------------------


void RendererSetLightCB(double r_light, double g_light, double b_light, 
                   double view_x, double view_y, double view_z)
/********************************************************************************/
/* Callback for modifications via Open Inventor light editor			*/
/********************************************************************************/
{
  int VP_LIGHT;
  switch (current_light){
    case 1 : VP_LIGHT = VP_LIGHT0;
         break;
    case 2 : VP_LIGHT = VP_LIGHT1;
         break;
    case 3 : VP_LIGHT = VP_LIGHT2;
         break;
    case 4 : VP_LIGHT = VP_LIGHT3;
         break;
    case 5 : VP_LIGHT = VP_LIGHT4;
         break;
    case 6 : VP_LIGHT = VP_LIGHT5;
         break;
    default: fprintf(stderr, "light number %d is not valid",current_light);
         exit(1);  
  };
  vpSetLight(vpc, VP_LIGHT, VP_DIRECTION, view_x, view_y,view_z);
  vpSetLight(vpc, VP_LIGHT, VP_COLOR, r_light, g_light, b_light);
  //vpSetLight(vpc, VP_LIGHT, VP_DIRECTION, x_direction, 
    //        y_direction, z_direction);
  vpShadeTable(vpc);
  multi_rotate_image(1,0.0,1);
}
// ---------------------------------------------------------------

void modify_light_numbers()
/********************************************************************************/
/* popups the menu for modifying the number of the light that can be modified 	*/
/* in modify_light_properties and for disabling/enabling the light		*/
/********************************************************************************/
{
  int light_enabled[6];
  vpGeti(vpc, VP_LIGHT0, &light_enabled[0]); 
  vpGeti(vpc, VP_LIGHT1, &light_enabled[1]); 
  vpGeti(vpc, VP_LIGHT2, &light_enabled[2]); 
  vpGeti(vpc, VP_LIGHT3, &light_enabled[3]); 
  vpGeti(vpc, VP_LIGHT4, &light_enabled[4]); 
  vpGeti(vpc, VP_LIGHT5, &light_enabled[5]); 
  light_numbers_menu(&current_light, &light_enabled[0]);
  //vpEnable(vpc, VP_LIGHT0, light_enabled[0]);
  //vpEnable(vpc, VP_LIGHT1, light_enabled[1]);
  //vpEnable(vpc, VP_LIGHT2, light_enabled[2]);
  //vpEnable(vpc, VP_LIGHT3, light_enabled[3]);
  //vpEnable(vpc, VP_LIGHT4, light_enabled[4]);
  //vpEnable(vpc, VP_LIGHT5, light_enabled[5]);
}
// ---------------------------------------------------------------

void setLights(int *light_enabled)
{
  vpEnable(vpc, VP_LIGHT0, light_enabled[0]);
  vpEnable(vpc, VP_LIGHT1, light_enabled[1]);
  vpEnable(vpc, VP_LIGHT2, light_enabled[2]);
  vpEnable(vpc, VP_LIGHT3, light_enabled[3]);
  vpEnable(vpc, VP_LIGHT4, light_enabled[4]);
  vpEnable(vpc, VP_LIGHT5, light_enabled[5]);
}

void activate_rotation_slider()
/************************************************************************/
/* popups the menu for modifying the rotation				*/
/************************************************************************/
{
  rotation_slider(&number_of_rotations, &rotation_increment);
}

// ---------------------------------------------------------------

void activate_histogram_slider()
/************************************************************************/
/* popups the menu for modifying the histogram boundaries		*/
/************************************************************************/
{
  histogram_slider(&histogram_minimum, &histogram_maximum);
}

// ---------------------------------------------------------------

void activate_scaling_sliders()
/************************************************************************/
/* popups the menu for modifying the rotation				*/
/************************************************************************/
{
  scale_voxels(&scale_x, &scale_y, &scale_z);
  printf("Scale by %f %f %f\n ", scale_x, scale_y, scale_z);
  multi_rotate_image(1,0.0,1);
}
// ---------------------------------------------------------------

void activate_depth_cueing_slider()
/************************************************************************/
/* popups the menu for modifying the rotation				*/
/************************************************************************/
{
  /* initial depth cueing values */ 
  static double front_factor=1.0, fog_density=2.0; 
  int enabled;
  /* is depth cueing enabled? */ 
  vpGeti(vpc, VP_DEPTH_CUE, &enabled); 
  depth_cueing_slider(&front_factor, &fog_density, (char)enabled);
  vpSetDepthCueing(vpc, front_factor, fog_density);
}
// ---------------------------------------------------------------

void setDepthCueingValues(double front_factor, double fog_density)
/************************************************************************/
/* sets new depth cueing values						*/
/************************************************************************/
{
  int enabled=0;
  vpGeti(vpc, VP_DEPTH_CUE, &enabled); 
  vpSetDepthCueing(vpc, front_factor, fog_density);
  if (!enabled) {
    show_message("Depth cueing is currently switched off",
                 "rendering is without depth cueing","");
  }
  multi_rotate_image(1,0.0,1);    
}

// ---------------------------------------------------------------

void setMaximumOpacity(double max_opacity)
/************************************************************************/
/* sets the maximum opacity (the threshold for early ray termination)	*/
/************************************************************************/
{
  vpSetd(vpc, VP_MAX_RAY_OPACITY, max_opacity);
}

// ---------------------------------------------------------------

void setMinimumOpacity(double min_opacity)
/************************************************************************/
/* sets the mimimum opacity for voxels					*/
/* below this threshold the voxels are transparent			*/
/************************************************************************/
{
  vpSetd(vpc, VP_MIN_VOXEL_OPACITY, min_opacity);
}

// ---------------------------------------------------------------

void activate_maximum_opacity_slider()
/************************************************************************/
/* popups the menu for modifying the maximum opacity			*/
/************************************************************************/
{
  double max_opacity;
  vpGetd(vpc, VP_MAX_RAY_OPACITY, &max_opacity);
  maximum_opacity_slider(&max_opacity);
  //vpSetd(vpc, VP_MAX_RAY_OPACITY, max_opacity);
  //multi_rotate_image(1, 0.0, 1);
}
// ---------------------------------------------------------------

void activate_minimum_opacity_slider()
/************************************************************************/
/* popups the menu for modifying the minimum opacity			*/
/************************************************************************/
{
  double min_opacity;
  vpGetd(vpc, VP_MIN_VOXEL_OPACITY, &min_opacity);
  minimum_opacity_slider(&min_opacity);
  vpSetd(vpc, VP_MIN_VOXEL_OPACITY, min_opacity);
  multi_rotate_image(1, 0.0, 1);
}

// ---------------------------------------------------------------

void calculate_histogram(){
/************************************************************************/
/* calculates the histogram of a volume file     			*/
/************************************************************************/
  int k, density_fd;
  unsigned density_size = data_xlen * data_ylen * data_zlen; 
  unsigned char *density, *tmp;
  unsigned int histogram_array[256];
  // load the density data 
  density = new unsigned char[density_size];
  if ((density_fd = open(density_filename, 0)) < 0) {
    show_message("Could not open scalar data file:",density_filename,"for histogram calculations");
    return;
  }
  if (lseek(density_fd, data_header, 0) < 0) {
    show_message("could not skip header of file:",density_filename,
    "for histogram calculations");
    return;
  }
  if (read(density_fd, density, density_size) != density_size) {
    show_message("could not read data form file:", density_filename,
                 "for histogram calculations");
    return;
  }
  // initialize histogram_array
  for (k =0; k < 256; k++){
    histogram_array[k] = 0;
  }
  // scan_buffer to create histogram
  tmp = density;
  for (k=0; k < density_size; k++){
    histogram_array[*density]++; density++;
  }
  // free buffer and close file
  density = tmp;
  delete density;
  close(density_fd);
  // draw histogram 
  draw_histogram(histogram_array, histogram_minimum, histogram_maximum);
}

// ---------------------------------------------------------------

void activate_xyplot(){
/************************************************************************/
/* calls popup menu for a new classification				*/
/************************************************************************/
  static unsigned char first_time = 1;
  if (first_time){ 
    int i;
    for (i = 0; i < n_scalar; i++){
      SRX[i]=(float)DensityRampX[i];
      SRY[i]=DensityRampY[i];
    }
    for (i = 0; i < n_gradient; i++){
      GRX[i]=(float)GradientRampX[i];
      GRY[i]=GradientRampY[i];
    }
    active_xyplot();
    first_time --;
  }
  else {
    active_xyplot();
  }
}
// ---------------------------------------------------------------

void new_classification(){
/************************************************************************/
/* do the new classification						*/
/************************************************************************/
    static int raw_volume_loaded = 0, octree_loaded = 0;
    int *ScalarRampX, *GradientRampX;
    int density_fd;	/* file descriptor for raw volume data (input) */
    int volume_fd, octree_fd, output_fd;	/* file descriptor for classified volume (output) */
    unsigned char *density; /* buffer for density data */
    unsigned density_size;  /* size of density data */
    int i;
    printf("n_scalar and n_gradient: %d %d\n",n_scalar,n_gradient);
    FORCEDRAW = 1;
    ScalarRampX = new int [n_scalar + 1];
    GradientRampX = new int [n_gradient + 1];
    for (i = 0; i < n_scalar; i++){
      ScalarRampX[i] = (int)(SRX[i]+0.49);
    }
    for (i = 0; i < n_gradient; i++){
      GradientRampX[i] = (int)(GRX[i]+0.49);
    }
    if (use_original_data) {
	/* allocate space for the raw data */
        printf("data lengths are %d %d %d\n",data_xlen,data_ylen,data_zlen);
        printf("data header length is %d \n",data_header);
	density_size = data_xlen * data_ylen * data_zlen;
        density = new unsigned char[density_size];
	if (density == NULL) {
          show_message("Exception during new classification:","out of memory","");
	  return;
	}

	/* load the raw data */
	if ((density_fd = open(density_filename, 0)) < 0) {
          show_message("could not open file:",density_filename,
                       "no new classification!");
	  return;
	}
	if (lseek(density_fd, data_header, 0) < 0) {
          show_message("could not read data from",density_filename,
                       "no new classification!");
	  return;
	}
	if (read(density_fd, density, density_size) != density_size) {
          show_message("could not read data from",density_filename,
                       "no new classification!");
          return;
	}
	close(density_fd);
    } else {
        if (!raw_volume_loaded){        
	/* load the unclassified volume data (the one with normals, etc.) */
        printf("loading raw_volume \n");
	if ((volume_fd = open(raw_volume_filename, 0)) < 0) {
          show_message("could not open file:",raw_volume_filename,
                       "no new classification!");
          return;
	}
	if (vpLoadRawVolume(vpc, volume_fd) != VP_OK) {
          show_message("VolPack error:",vpGetErrorString(vpGetError(vpc)),
            "no new classification");
	  return;
	}
        raw_volume_loaded = 1;
	close(volume_fd);
        }
    }
    vpSetVolumeSize(vpc, data_xlen, data_ylen, data_zlen);
    vpSetVoxelSize(vpc, BYTES_PER_VOXEL, VOXEL_FIELDS, SHADE_FIELDS, CLSFY_FIELDS);
    vpSetVoxelField(vpc, NORMAL_FIELD, NORMAL_SIZE, NORMAL_OFFSET, NORMAL_MAX);
    vpSetVoxelField(vpc, DENSITY_FIELD, DENSITY_SIZE, DENSITY_OFFSET, DENSITY_MAX);
    vpSetVoxelField(vpc, GRADIENT_FIELD, GRADIENT_SIZE, GRADIENT_OFFSET, GRADIENT_MAX);

    /* set the classification function */
    vpRamp(density_ramp, sizeof(float), n_scalar, ScalarRampX,
	       SRY);
    vpSetClassifierTable(vpc, DENSITY_PARAM, DENSITY_FIELD, density_ramp,
			     sizeof(density_ramp));
    vpRamp(gradient_ramp, sizeof(float), n_gradient,
	       GradientRampX, GRY);
    vpSetClassifierTable(vpc, GRADIENT_PARAM, GRADIENT_FIELD,
			     gradient_ramp, sizeof(gradient_ramp));
    /* load the octree */
    if (use_octree && (!octree_loaded)) {
      /* load the octree */
      printf("Loading the octree!\n");
      if ((octree_fd = open(octree_filename, 0)) < 0) {
        show_message("could not open octree file:",octree_filename,
                       "no new classification!");
        return;
      }
      if (vpLoadMinMaxOctree(vpc, octree_fd) != VP_OK) {
        show_message("VolPack error:",vpGetErrorString(vpGetError(vpc)),
                       "no new classification!");
        return;
      }
      octree_loaded = 1;
      close(octree_fd);
    }
    /* classify */
    if (use_original_data) {
	if (vpClassifyScalars(vpc, density, density_size, DENSITY_FIELD,
			      GRADIENT_FIELD, NORMAL_FIELD) != VP_OK) {
          show_message("VolPack error:",vpGetErrorString(vpGetError(vpc)),
          "no new classification");
          return;
	}
    } else {
	if (vpClassifyVolume(vpc) != VP_OK) {
	    show_message("VolPack error",
		    vpGetErrorString(vpGetError(vpc)),
                    "no new classifcation");
	    return;
	}
    }
    /* store the classified volume */
    /* FG - do not store, if there are two dots in current_filename */
    if (strstr( strstr(current_filename,".") + 1, "." )==NULL) {
      if ((output_fd = creat(current_filename, 0644)) < 0) {
        show_message("could not open",current_filename,
                     "as a new classified volume file");
        return;
      }
      if (vpStoreClassifiedVolume(vpc, output_fd) != VP_OK) {
	  show_message("VolPack error:",
		vpGetErrorString(vpGetError(vpc)),"during new classification");
	  exit(1);
      }
      close(output_fd);
    }

    /* free memory */
    if (use_original_data) delete density;
    delete ScalarRampX;  delete GradientRampX;
    /* load new classified volume */
    /* FG - do reload, if there are two dots in current_filename */
    if (strstr( strstr(current_filename,".") + 1, "." )==NULL) {
      load_classified_volume(current_filename);
    }
    /* render volume */
    multi_rotate_image(1, 0.0, 1);
}
// ---------------------------------------------------------------

void StorePGM(char *output_filename)
/************************************************************************/
/* Store the image in Josefs Pankofers PGM or PPM format		*/
/************************************************************************/
{
    FILE *image_fp;	/* file descriptor for image (output) */

/* P5 is for luminance only, P6 for rgb storage */
#define PGM_MAGIC1	'P'
#define RPGM_MAGIC5	'5'
#define RPGM_MAGIC6	'6'
  if ((image_fp = fopen(output_filename, "w")) == NULL) {
    show_message("could not open output file",output_filename,"");
    return;
  }
  switch (CurrentColorMode){
    case LUMINANCE_MODE: 
      fprintf(image_fp, "%c%c\n%d %d\n%d\n", PGM_MAGIC1, RPGM_MAGIC5,
	        image_width, image_height, 255);
      fwrite(image[0], 1, image_width*image_height, image_fp);
      break;
    case RGB_MODE:
      fprintf(image_fp, "%c%c\n%d %d\n%d\n", PGM_MAGIC1, RPGM_MAGIC6,
	        image_width, image_height, 255);
      fwrite(image[0], 1, image_width*image_height*3, image_fp);
      break;
    default: 
      fprintf(stderr,"StorePGM Error: Unknown CurrentColorMode");
      exit(1);
  }
  fclose(image_fp);
}
// ---------------------------------------------------------------

unsigned char *ImgRenderCB(const ImageInfo &ifo, int PicNo)
/************************************************************************/
/* Render Callback for ImgViewer Callbacks	       			*/
/* returns the image adress 						*/
/************************************************************************/
{
  vpMatrix4 render_matrix;
  //int image_width, image_height;
  static int called_first_time = 1;
  static float old_eye_distance, current_eye_distance; 
  static float zoom_factor;
  float camera_zoom_factor=1.0; 
  vpCurrentMatrix(vpc, VP_MODEL);
  vpIdentityMatrix(vpc);
  //vpCurrentMatrix(vpc, VP_VIEW);
  vpGetMatrix(vpc, VP_MODEL, render_matrix);
  
  image_width = ifo.SizeX;
  image_height = ifo.SizeY;
  switch (CurrentColorMode){
    case LUMINANCE_MODE: vpSetImage(vpc, image[PicNo], image_width, image_height, 
                                           image_width, VP_LUMINANCE);
                      break;
    case RGB_MODE: vpSetImage(vpc, image[PicNo], image_width, image_height, 
                                           image_width*3, VP_RGB);
                      break;
    default: fprintf(stderr,"Error: Unknown CurrentColorMode");
                      exit(1);
  }
  /* seems no one is using the same convention for axis ! */
  //printf("4 * 4 matrix is:\n");
  //printf("%f %f %f %f\n",ifo.perspect[0],ifo.perspect[1],ifo.perspect[2],ifo.perspect[3]);
  //printf("%f %f %f %f\n",ifo.perspect[4],ifo.perspect[5],ifo.perspect[6],ifo.perspect[7]);
  //printf("%f %f %f %f\n",ifo.perspect[8],ifo.perspect[9],ifo.perspect[10],ifo.perspect[11]);
  //printf("%f %f %f %f\n",ifo.perspect[12],ifo.perspect[13],ifo.perspect[14],ifo.perspect[15]);
  render_matrix[0][0] = ifo.perspect[0];
  render_matrix[1][0] = ifo.perspect[1];
  render_matrix[2][0] = -ifo.perspect[2];
  render_matrix[0][3] = ifo.perspect[3];

  render_matrix[0][1] = ifo.perspect[4];
  render_matrix[1][1] = ifo.perspect[5];
  render_matrix[2][1] = -ifo.perspect[6];
  render_matrix[1][3] = ifo.perspect[7];

  render_matrix[0][2] = -ifo.perspect[8];
  render_matrix[1][2] = -ifo.perspect[9];
  render_matrix[2][2] = ifo.perspect[10];
  render_matrix[2][3] = -ifo.perspect[11];

  render_matrix[3][0] = ifo.perspect[12];
  render_matrix[3][1] = ifo.perspect[13];
  render_matrix[3][2] = ifo.perspect[14];
  render_matrix[3][3] = ifo.perspect[15];

  render_matrix[0][3] = -ifo.perspect[12]/ifo.perspect[14];
  render_matrix[1][3] = -ifo.perspect[13]/ifo.perspect[14];
  render_matrix[3][0] = 0.;
  render_matrix[3][1] = 0.;
  render_matrix[3][2] = 0.;
  render_matrix[3][3] = 1.;
  //printf("4 * 4 volpack matrix is:\n");
  //printf("%f %f %f %f\n",render_matrix[0][0],render_matrix[0][1],
  //                       render_matrix[0][2],render_matrix[0][3]);
  // printf("%f %f %f %f\n",render_matrix[1][0],render_matrix[1][1],
  //                       render_matrix[1][2],render_matrix[1][3]);
  //printf("%f %f %f %f\n",render_matrix[2][0],render_matrix[2][1],
  //                       render_matrix[2][2],render_matrix[2][3]);
  //printf("%f %f %f %f\n",render_matrix[3][0],render_matrix[3][1],
  //                       render_matrix[3][2],render_matrix[3][3]);


  vpSetMatrix(vpc, render_matrix);

  if (called_first_time) {
    called_first_time = 0;
    old_eye_distance = ifo.eye[2]*ifo.eye[2] +
                       ifo.eye[1]*ifo.eye[1] + 
                       ifo.eye[0]*ifo.eye[0]; 
    zoom_factor = 1;
  }
  else {
    current_eye_distance = ifo.eye[2]*ifo.eye[2] +
                       ifo.eye[1]*ifo.eye[1] + 
                       ifo.eye[0]*ifo.eye[0]; 
    zoom_factor = old_eye_distance/current_eye_distance*zoom_factor;
// this is a preliminary solution for the viewer makes a difference
// between orthographic and perspective camera, until now
// in volpack there is only orthographic projection
// global.Camera_Type = 1 : Perspective Camera
// global.Camera_Type = 0 : Orthographic Camera
    if (global.Camera_Type){ 
      camera_zoom_factor = global.viewer->getCameraZoom()/90.0;
    }
    else{
      camera_zoom_factor = global.viewer->getCameraZoom()/6.228001;
    }
    vpSeti(vpc, VP_CONCAT_MODE, VP_CONCAT_RIGHT);
    vpScale(vpc, (double)zoom_factor*scale_x/camera_zoom_factor, 
                 (double)zoom_factor*scale_y/camera_zoom_factor, 
                 (double)zoom_factor*scale_z/camera_zoom_factor);  
    vpSeti(vpc, VP_CONCAT_MODE, VP_CONCAT_LEFT);
    old_eye_distance = current_eye_distance;
  }

  switch(global.viewer->getDrawStyle(ImgViewer::STILL)){
    case 0:
    // volpack shading
      vpc->shading_mode = LOOKUP_SHADER;
      if (vpShadeTable(vpc) != VP_OK) {
        fprintf(stderr, "VolPack error calling vpShadeTable: %s\n",
        vpGetErrorString(vpGetError(vpc)));
        exit(1);
      }
      break;
    case 1:
    // callback shading
      if (CurrentColorMode == RGB_MODE){
        show_message("The implemented shader callback works only if",
                     "RGB colors are not used for rendering!",
                     "Continue using the VolPack Lookup Shader");
        vpc->shading_mode = LOOKUP_SHADER;
        if (vpShadeTable(vpc) != VP_OK) {
          fprintf(stderr, "VolPack error calling vpShadeTable: %s\n",
          vpGetErrorString(vpGetError(vpc)));
          exit(1);
        }
        break;
      };
      vpc->shade_func = (void (*)())callback_shader;
      vpc->shading_mode = CALLBACK_SHADER;
      break;
  }

  /* call the renderer only if the DrawFlag is not zero */
  if (DrawFlag){
    if (use_clvolume) {
	if (vpRenderClassifiedVolume(vpc) != VP_OK) {
	    fprintf(stderr, "VolPack vpRenderClassifiedVolume error: %s\n",
		vpGetErrorString(vpGetError(vpc)));
	    //exit(1);
	}
     } else {
	if (vpRenderRawVolume(vpc) != VP_OK) {
	    fprintf(stderr, "VolPack error vpRenderRawVolume: %s\n",
		    vpGetErrorString(vpGetError(vpc)));
	    exit(1);
	}
    }
  }
  return (unsigned char *)image[PicNo];
}

// ---------------------------------------------------------------
void store_materials(const char *material_filename){
  double r_ambient, r_diffuse, r_specular,
         g_ambient, g_diffuse, g_specular,
         b_ambient, b_diffuse, b_specular,
         shinyness, dummy_shinyness;
  int i,k;
  FILE *material_fp;
  if ((material_fp = fopen(material_filename, "w")) == NULL) {
    show_message("could not open material file:",material_filename,
    "for storage of material characteristics");
    return;
  }
  switch (CurrentColorMode){
    case LUMINANCE_MODE: fprintf(material_fp,"0\n");
                         break;
    case RGB_MODE: fprintf(material_fp,"1\n");
                   break;
    default: fprintf(stderr,"Error: Unknown CurrentColorMode\n");
                      exit(1);
  }
  fprintf(material_fp,"%d\n",NUM_MATERIALS);
  // foreach active material
  for(i=0;i<NUM_MATERIALS;i++){
    fprintf(material_fp,"%d\n",material_weight_points[i]);
    // store the material weight points
    for(k=0;k<material_weight_points[i];k++){
      fprintf(material_fp,"%f %f\n",WeightPointsX[i][k],WeightPointsY[i][k]); 
    }
    // get the color parameters
    switch (i){
      case 0 : VP_MATERIAL = VP_MATERIAL0;
           break;
      case 1 : VP_MATERIAL = VP_MATERIAL1;
           break;
      case 2 : VP_MATERIAL = VP_MATERIAL2;
           break;
      case 3 : VP_MATERIAL = VP_MATERIAL3;
           break;
      case 4 : VP_MATERIAL = VP_MATERIAL4;
           break;
      case 5 : VP_MATERIAL = VP_MATERIAL5;
           break;
      default: fprintf(stderr, "store_materials: material number %d is not valid\n",i);
           exit(1);  
    };
    vpGetMaterial(vpc, VP_MATERIAL, VP_AMBIENT, VP_EXTERIOR, 
                   &r_ambient, &g_ambient, &b_ambient);
    vpGetMaterial(vpc, VP_MATERIAL, VP_DIFFUSE, VP_EXTERIOR, 
                   &r_diffuse, &g_diffuse, &b_diffuse);
    vpGetMaterial(vpc, VP_MATERIAL, VP_SPECULAR, VP_EXTERIOR, 
                   &r_specular, &g_specular, &b_specular);
    vpGetMaterial(vpc, VP_MATERIAL, VP_SHINYNESS, VP_EXTERIOR, 
                   &shinyness, &dummy_shinyness, &dummy_shinyness);
    // store the color parameters
      fprintf(material_fp,"%f %f %f\n", (float)r_ambient, 
              (float)g_ambient, (float)b_ambient);
      fprintf(material_fp,"%f %f %f\n", (float)r_diffuse, 
              (float)g_diffuse, (float)b_diffuse);
      fprintf(material_fp,"%f %f %f\n", (float)r_ambient, 
              (float)g_specular, (float)b_specular);
      fprintf(material_fp,"%f\n",(float)shinyness);
  }
  fclose(material_fp);
}
// ---------------------------------------------------------------



void load_materials(const char *material_filename){
  double r_ambient, r_diffuse, r_specular,
         g_ambient, g_diffuse, g_specular,
         b_ambient, b_diffuse, b_specular,
         shinyness;
  int i,k, dummy_color_mode;
  FILE *material_fp;
  if ((material_fp = fopen(material_filename, "r")) == NULL) {
    show_message("could not open material file:", material_filename,"");
    return;
  }
  fscanf(material_fp,"%d",&dummy_color_mode);
  if ((dummy_color_mode != 0) && (dummy_color_mode != 1)){
    printf("dummy_color_mode us %d\n",dummy_color_mode);
    show_message("May be an older file version","in the first line  of the file",
                 "must be either 0 or 1");
    return;
  }
  fscanf(material_fp,"%d",&NUM_MATERIALS);
  if (NUM_MATERIALS > MAX_MATERIALS){
    show_message("This number of materials","is not supported!","");
    fclose(material_fp);
    return;
  }
  update_material_numbers(NUM_MATERIALS);
  // foreach active material
  for(i=0;i<NUM_MATERIALS;i++){
    fscanf(material_fp,"%d",&(material_weight_points[i]));
    // store the material weight points
    for(k=0;k<material_weight_points[i];k++){
      fscanf(material_fp,"%f %f",&(WeightPointsX[i][k]),&(WeightPointsY[i][k])); 
      printf("material points are %f %f\n",WeightPointsX[i][k],WeightPointsY[i][k]); 
    }
    // get the color parameters
      fscanf(material_fp,"%lf %lf %lf", &r_ambient, 
              &g_ambient, &b_ambient);
      printf("r_ambient,b_ambient,g_ambient %lf %lf %lf\n",r_ambient,b_ambient,g_ambient);
      fscanf(material_fp,"%lf %lf %lf", &r_diffuse, 
              &g_diffuse, &b_diffuse);
      printf("r_diffuse,b_diffuse,g_diffuse %lf %lf %lf\n",r_diffuse,b_diffuse,
              g_diffuse);
      fscanf(material_fp,"%lf %lf %lf", &r_specular, 
              &g_specular, &b_specular);
      fscanf(material_fp,"%lf",&shinyness);
    // set the color parameters
    switch (i){
      case 0 : VP_MATERIAL = VP_MATERIAL0;
           break;
      case 1 : VP_MATERIAL = VP_MATERIAL1;
           break;
      case 2 : VP_MATERIAL = VP_MATERIAL2;
           break;
      case 3 : VP_MATERIAL = VP_MATERIAL3;
           break;
      case 4 : VP_MATERIAL = VP_MATERIAL4;
           break;
      case 5 : VP_MATERIAL = VP_MATERIAL5;
           break;
      default: fprintf(stderr, "load_materials: NUM_MATERIALS %d is not valid\n",i);
           exit(1);  
    };
    vpSetMaterial(vpc, VP_MATERIAL, VP_AMBIENT, VP_BOTH_SIDES, 
                   r_ambient, g_ambient, b_ambient);
    vpSetMaterial(vpc, VP_MATERIAL, VP_DIFFUSE, VP_BOTH_SIDES, 
                   r_diffuse, g_diffuse, b_diffuse);
    vpSetMaterial(vpc, VP_MATERIAL, VP_SPECULAR, VP_BOTH_SIDES, 
                   r_specular, g_specular, b_specular);
    vpSetMaterial(vpc, VP_MATERIAL, VP_SHINYNESS, VP_BOTH_SIDES, 
                   shinyness, shinyness, shinyness);
  }
  fclose(material_fp);
  set_material_ramps();
}
