#ifndef _VIEWER_H
#define _VIEWER_H

void initViewer(ImgExaminerViewer* &, Widget, Widget);
void ModifyLightCB(double r_light, double g_light, double b_light,
              double x_direction, double y_direction, double z_direction);
void RendererSetLightCB(double r_light, double g_light, double b_light, 
                   double x_direction, double y_direction, double z_direction);

#endif
