#ifndef _RENDER_H
#define _RENDER_H

void initialize_volpack();
void multi_rotate_image(int axis, float angle, int number_of_rotations);
unsigned char *ImgRenderCB(const ImageInfo &ifo, int soc);
unsigned char *rotate_image(int axis, float angle);
void set_image_sequence_name(char *sequence_filename);
void load_and_render_sequence(char *sequence_filename);
void load_unclassified_sequence(char *sequence_filename);
void load_classified_volume(char *clVOLUME_FILE);
void load_and_classify_volume(char *filename);
void load_and_classify_next_volume(char *filename);
void load_image(char *filename);
void change_filenames(char *filename);
void switch_image_store_mode(void);
void switch_color_mode(void);
void switch_octree_mode(void);
void build_octree(void);
void build_raw_volume(void);
void vpToggle(int vp_option);
void modify_material_properties(void);
void modify_material_properties_applyCB(double r_ambient, double g_ambient, double b_ambient,
		       double r_diffuse,double g_diffuse,double b_diffuse,
                       double r_specular,double  g_specular,double b_specular,
		       double shinyness);
void update_material_properties(
  double *r_ambient, double *g_ambient, double *b_ambient,
  double *r_diffuse, double *g_diffuse, double *b_diffuse,
  double *r_specular, double *g_specular, double *b_specular,
  double *shinyness, int new_current_material, int *material_weight_points);
void update_material_numbers(int number_of_materials);
void modify_material_numbers(void);
void setLights(int *light_enabled);
void setDepthCueingValues(double front_factor, double fog_density);
void setMaximumOpacity(double max_opacity);
void setMinimumOpacity(double min_opacity);
void setLightDirection(double x_direction, double y_direction, double z_direction);
void setLightColors(double r_light, double g_light, double b_light);
void set_material_ramp(int material_weight_points);
void set_material_ramps();
void modify_light_properties(void);
void modify_light_numbers(void);
void activate_histogram_slider(void);
void activate_rotation_slider(void);
void activate_scaling_sliders(void);
void activate_maximum_opacity_slider(void);
void activate_minimum_opacity_slider(void);
void activate_depth_cueing_slider(void);
void activate_xyplot(void);
void new_classification(void);
void StorePGM(char *filename);
void RendererSetLightCB(double r_light, double g_light, double b_light, 
                   double x_direction, double y_direction, double z_direction);
void switch_light_editor(void);
void calculate_histogram();
void store_materials(const char *material_filename);
void load_materials(const char *material_filename);
int get_number_of_materials();
void set_material_weight_points(int new_material_weight_points,int current_material);
#endif
