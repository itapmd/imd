void calculate_histogram(){
histogram_array[256];
int k, density_fd;
unsigned density_size = data_xlen * data_ylen * data_zlen;
unsigned char *density = new unsigned char [density_size];
unsigned int histogramm_array[256];
/* load the density data */
if ((density_fd = open(density_filename, 0)) < 0) {
    perror("open");
    fprintf(stderr, "could not open %s\n", density_filename);
    exit(1);
}
if (lseek(density_fd, data_header, 0) < 0) {
    perror("seek");
    fprintf(stderr, "could not skip header of file %s\n", density_filename);
    exit(1);
}
if (read(density_fd, density, density_size) != density_size) {
    perror("read");
    fprintf(stderr, "could not read data from %s\n", density_filename);
    exit(1);
}
// initialize histogram_array
for (k =0; k < 256; k++){
  histogram_array[k] = 0;
}
// scan_buffer to create histogram
for (k=0; k < density; k++){
  histogramm_array[*density]++; density++;
}
// free buffer and close file
delete density;
close(density_fd);
// draw histogram 
draw_histogram(histogram_array);
}

void draw_histogram(int* histogram_array){
// simple xyplot
  typedef struct{
    FL_FORM *histogram_form;
    FL_OBJECT *histogram_xyplot;
    FL_OBJECT *histogram_quit;  
  } FD_histogram;
  float scalar_value[256], counts[256];
  int k;

  FL_OBJECT *obj;
  FD_histogram *fdui = (FD_histogram *) fl_calloc(1, sizeof(FD_histogram));
  if (first_time) forms_initialization();
  for (k=0;k<256;k++){
    scalar_value[k] = (float) k;
    counts[k] = (float) histogram_array[k];
  }
  fdui->histogram_form = fl_bgn_form(FL_NO_BOX, 560, 310);
    obj = fl_add_box(FL_UP_BOX,0,0,560,310,"");
    fdui->histogram_xyplot = obj = fl_add_xyplot(FL_IMPULSE_XYPLOT,40,90,460,220,"Histogram");
    fl_set_object_boxtype(obj,FL_UP_BOX);
    fdui->histogram_quit = obj = fl_add_button(FL_NORMAL_BUTTON,430,20,70,30,"Quit");
  fl_end_form();

  fl_set_xyplot_data(fdui->histogram_xyplot, scalar_value, counts, 256, "","","");
  fl_set_xyplot_xtics(fdui->histogram_xyplot, 14, 3);
  fl_set_xyplot_ytics(fdui->histogram_xyplot, 14, 4);

  /* show the form */
  fl_show_form(fdui->histogram_xyplot,FL_PLACE_MOUSE,FL_TRANSIENT,"Histogram");
  while (obj != fdui->histogram_quit){
    obj = fl_do_forms();
  }
  fl_hide_form(fdui->histogram_xyplot);
  fl_free_form(fdui->histogram_xyplot);
}
