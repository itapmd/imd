// **************************************************************************
//
// Description:   IMD plugin for simulation steering with IMD (ITAP)
//
// Author:        Juergen Schulze
//
// Creation Date: 2003-01-09
//
// **************************************************************************

#include <coVRModuleSupport.h>

#include <Performer/pf/pfGeode.h>
#include <Performer/pf/pfDCS.h>
#include <Performer/pr/pfMaterial.h>

#include "AtomScenegraph.h"

/// Constructor
AtomScenegraph::AtomScenegraph()
{
  const char* trianglesString;

  material = NULL;
  geoState = NULL;
  numAtoms = 0;

  trianglesString = CoviseConfig::getEntry("IMDPluginConfig.AtomTriangles");
  if (trianglesString)
  {
    sphereTriangles = atoi(trianglesString);
    if (sphereTriangles<1) sphereTriangles = 1;
  }
  else
  {
    sphereTriangles = 20;
  }

  createMaterial();
  createGeoState();
  sceneDCS = new pfDCS();
  
  covise = cover->getObjectsRoot();
  covise->addChild(sceneDCS);
}

/// Destructor
AtomScenegraph::~AtomScenegraph()
{
  covise->removeChild(sceneDCS);
  pfDelete(sceneDCS);
}

/** Add a sphere (atom) to the scene.
  @param xpos,ypos,zpos  location of atom center [world coordinates]
  @param size            diameter of the atom [mm]
  @param color           atom color on color gradient [0..1]
*/
void AtomScenegraph::addAtom(float xpos, float ypos, float zpos, 
  float size, float red, float green, float blue)
{
  pfGeode*  sphereGeode;
  pfGeoSet* sphereGeoSet;
  pfVec4    *sphereColor;
  pfMatrix  sphereMatrix;

  // Create nodes:
  sphereGeode = new pfGeode();

  // Create sphere:
  sphereGeoSet = pfdNewSphere(sphereTriangles, pfGetSharedArena());
  sphereGeoSet->setGState(geoState);
  
  // Set color:
  sphereColor = new pfVec4();
  sphereColor->set(red, green, blue, 1.0);
  sphereGeoSet->setAttr(PFGS_COLOR4, PFGS_OVERALL, sphereColor, NULL);
  
  // Set size:
  sphereMatrix.makeScale(size, size, size);
  pfdXformGSet(sphereGeoSet, sphereMatrix);

  // Set position:
  sphereMatrix.makeTrans(xpos, ypos, zpos);
  pfdXformGSet(sphereGeoSet, sphereMatrix);

  // Finalize:
  sphereGeode->addGSet(sphereGeoSet);
  sceneDCS->addChild(sphereGeode);
  
  ++numAtoms;
}

/** Create the material which is to be used for all atoms;
*/
void AtomScenegraph::createMaterial()
{
  material = new pfMaterial();
  material->setSide(PFMTL_FRONT);
  material->setColorMode(PFMTL_FRONT, PFMTL_CMODE_AMBIENT_AND_DIFFUSE);
  material->setColor(PFMTL_AMBIENT,  0.2f, 0.2f, 0.2f);
  material->setColor(PFMTL_DIFFUSE,  0.9f, 0.9f, 0.9f);
  material->setColor(PFMTL_SPECULAR, 0.9f, 0.9f, 0.9f);
  material->setColor(PFMTL_EMISSION, 0.0f, 0.0f, 0.0f);
  material->setShininess(16.0f);
}

void AtomScenegraph::createGeoState()
{
  geoState = new pfGeoState();
  geoState->makeBasic();
  geoState->setAttr(PFSTATE_FRONTMTL, material);
  geoState->setMode(PFSTATE_ENLIGHTING, PF_ON);
  geoState->setMode(PFSTATE_CULLFACE, PFCF_BACK);
  geoState->setMode(PFSTATE_TRANSPARENCY, PFTR_OFF);
}

// EOF
