void set_dimensions(int *data_header, int *data_xlen,int *data_ylen,
			int *data_zlen, unsigned char *apply_button_pushed);
void depth_cueing_slider(double *front_factor, double *fog_density, char enabled);
void histogram_slider(int *histogram_minimum, int *histogram_maximum);
void rotation_slider(int *number_of_rotations, double *rotation_increment);
void show_message(const char *s1, const char *s2, const char *s3);
int show_question(const char *s1);
void maximum_opacity_slider(double *maximum_opacity);
void minimum_opacity_slider(double *minimum_opacity);
void scale_voxels(double *scale_x, double *scale_y, double *scale_z);
void active_xyplot(void);
void modify_lighting_popup(double *r_light, double *g_light, double *b_light,
                           double *x_direction, double *y_direction, 
                           double *z_direction, int current_light);
void modify_material_popup(double* r_ambient, double*g_ambient, double*b_ambient,
		       double*r_diffuse,double *g_diffuse,double *b_diffuse,
                       double*r_specular,double *g_specular,double *b_specular,
		       double*shinyness, 
                       int *current_material, int *number_of_materials,
                       int material_weight_points);
void material_numbers_sliders(int *material_count, int *current_material);
void light_numbers_menu(int *current_light, int *light_enabled);
void ask_server_name();
void ask_port_address();
void draw_histogram(unsigned int* histogram_array, 
                    int histogram_min, int histogram_max);
void remove_histogram();
void selectClassifiedFile(const char *filename);
void selectUnclassifiedFile(const char *filename);
void storeImage(const char *filename);
void selectStoreSequence(const char *filename);
void selectClassifiedLoadSequence(const char *filename);
void selectUnclassifiedLoadSequence(const char *filename);
