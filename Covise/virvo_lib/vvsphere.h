//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVSPHERE_H_
#define _VVSPHERE_H_

#include <stdlib.h>
#include "vvvecmath.h"

/** Vertex definition for vvSphere.
    @author Daniel Weiskopf
    @see vvSphere
*/
class vvVertex 
{
  public:
    vvVertex() {};

    vvVertex(const vvVertex&);
    const vvVertex & operator=(const vvVertex &);
    void  scale(float);

    float x;
    float y;
    float z;
    float th;
    float ph;
};

/** Texture coordinate for vvSphere.
  @author Daniel Weiskopf
  @see vvSphere
*/
class vvTexCoord 
{
  public:
    vvTexCoord() {};
    float t[3];
};

/** Triangle definition for vvSphere.
  @author Daniel Weiskopf
  @see vvSphere
*/
class vvTriangle 
{
  public:
    vvTriangle();
    vvTriangle(const vvTriangle&);
    const vvTriangle& operator=(const vvTriangle&);

    int v1;
    int v2;
    int v3;
    int visibility;
};

/** Generates texture coordinates for random spheres.
  @author Daniel Weiskopf
*/
class vvSphere
{
  public:
    vvSphere();     // Default constructor
    void render();
    void renderWireframe(int type=0);
    void initDodecaeder();
    void subdivide();
    void performCulling();
    void calculateTexCoords();
    void setRadius(float);
    void setVolumeDim(vvVector3*);
    void setViewMatrix(vvMatrix4*);
    void setTextureOffset(float*);

  private:
    static int dodecaederConnectivity[60];
    vvVertex*   mVertices;
    vvVertex*   mVerticesWorld;
    int         mNumVertices;
    vvTriangle* mTriangles;
    int         mNumTriangles;
    vvTexCoord* mTexCoords;
    double      mRadius;
    vvVector3    mDimCube;
    vvMatrix4    mModelView;
    float       texOffset[3];
    
    bool     isVisibleVertex(int vert);
    vvVertex midpoint(const vvVertex& a, const vvVertex& b) const;
    void     copyVerticesWorld();
};

#endif
