// **************************************************************************
//
// Description:   IMD plugin for simulation steering with IMD (ITAP)
//
// Author:        Juergen Schulze
//
// Creation Date: 2003-01-09
//
// **************************************************************************

#ifndef ATOM_SCENEGRAPH_H
#define ATOM_SCENEGRAPH_H

/** Encapsulation of Performer scenegraph related routines
  to create a scene of atoms.
*/
class AtomScenegraph
{
  private:
    pfDCS* sceneDCS;        ///< DCS for entire atom scene
    pfMaterial* material;   ///< Performer material for all atoms
    pfGeoState* geoState;   ///< Performer geo state for all atoms
    pfGroup* covise;        ///< add to covise to use the menus
    int sphereTriangles;    ///< number of triangles for spheres
    int numAtoms;           ///< number of atoms in scenegraph
    
    void createMaterial();
    void createGeoState();

  public:
    AtomScenegraph();
    ~AtomScenegraph();
    void addAtom(float, float, float, float, float, float, float);
};

#endif

// EOF
