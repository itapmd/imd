/**************************************************************************\
**                                                           (C)2002 RUS  **
**                                                                        **
** Description: Read IMD checkpoint files from ITAP.                      **
**                                                                        **
**                                                                        **
**                                                                        **
**                                                                        **
**                                                                        **
** Author:                                                                **
**                                                                        **
**                     Juergen Schulze-Doebold                            **
**     High Performance Computing Center University of Stuttgart          **
**                         Allmandring 30                                 **
**                         70550 Stuttgart                                **
**                                                                        **
** Cration Date: 03.09.2002                                               **
\**************************************************************************/

#include <coModule.h>
#include <vvtokenizer.h>
#include <vvtoolshed.h>
#include <vvarray.h>
#include <limits.h>
#include <float.h>
#include "ReadIMD.h"

/// Startup routine
int main(int argc, char *argv[])
{
  coReadIMD *application = new coReadIMD();
  application->start(argc, argv);
  return 0;
}

/// Constructor
coReadIMD::coReadIMD()
{
  set_module_description("Read IMD checkpoint files to create a list of points and scalar parameters.");

  // Create ports:
  poPoints = addOutPort("Location", "Set_Points", "Atom locations");
  poPoints->setInfo("Atom locations");

  poType = addOutPort("Type", "Set_Unstructured_S3D_Data", "Atom types");
  poType->setInfo("Atom types");

  poMass = addOutPort("Mass", "Set_Unstructured_S3D_Data", "Atom mass");
  poMass->setInfo("Atom mass");

  poSpeed = addOutPort("Speed", "Set_Unstructured_V3D_Data", "Atom speed");
  poSpeed->setInfo("Atom speed");

  poEpot = addOutPort("Epot", "Set_Unstructured_S3D_Data", "Epot");
  poEpot->setInfo("Epot");

  // Create parameters:
  pbrCheckpointFile = addBrowserParam("FilePath", "First IMD checkpoint file of sequence");
  pbrCheckpointFile->setValue("data/","*.cpt");

  pboPeriodic = addBooleanParam("PeriodicBoundaries", "Select to use periodic boundaries");
  pboPeriodic->setValue(false);

  pfvPeriodic = addFloatVectorParam("Splitter", "Splitter location for periodic boundaries");
  pfvPeriodic->setValue(0.0f, 0.0f, 0.0f);

  pboConstrainSpeed = addBooleanParam("ConstrainSpeed", "Constrain atom display with respect to speed (see MinSpeed and MaxSpeed)");
  pboConstrainSpeed->setValue(false);

  pfsSpeedMin = addFloatScalarParam("MinSpeed", "Minimum absolute speed for atoms to be displayed");
  pfsSpeedMin->setValue(0.0f);
  
  pfsSpeedMax = addFloatScalarParam("MaxSpeed", "Maximum absolute speed for atoms to be displayed");
  pfsSpeedMax->setValue(1.0f);

  pfsTypeOffset = addFloatScalarParam("TypeOffset", "Offset to add to type value");
  pfsTypeOffset->setValue(0.0f);

  pboIgnoreWarnings = addBooleanParam("IgnoreWarnings", "Ignore warnings when reading files");
  pboIgnoreWarnings->setValue(false);
}

/// @return absolute value of a vector
float coReadIMD::absVector(float x, float y, float z)
{
  return sqrt(x * x + y * y + z * z);
}

/// Compute routine: load checkpoint file
int coReadIMD::compute()
{
  const int DEF_DIM = 10000;
  vvArray<DO_Points*> aPoints(1,5);
  vvArray<DO_Unstructured_S3D_Data*> aType(1,5);
  vvArray<DO_Unstructured_S3D_Data*> aMass(1,5);
  vvArray<DO_Unstructured_V3D_Data*> aSpeed(1,5);
  vvArray<DO_Unstructured_S3D_Data*> aEpot(1,5);
  vvArray<float> number(0,DEF_DIM);
  vvArray<float> type(0,DEF_DIM);
  vvArray<float> mass(0,DEF_DIM);
  vvArray<float> epot(0,DEF_DIM);
  vvArray<float> x(0,DEF_DIM);
  vvArray<float> y(0,DEF_DIM);
  vvArray<float> z(0,DEF_DIM);
  vvArray<float> vx(0,DEF_DIM);
  vvArray<float> vy(0,DEF_DIM);
  vvArray<float> vz(0,DEF_DIM);
  DO_Points* doPoints = NULL;
  DO_Unstructured_S3D_Data* doType = NULL;
  DO_Unstructured_S3D_Data* doMass = NULL;
  DO_Unstructured_V3D_Data* doSpeed = NULL;
  DO_Unstructured_S3D_Data* doEpot = NULL;
  DO_Set* setPoints = NULL;
  DO_Set* setType = NULL;
  DO_Set* setMass = NULL;
  DO_Set* setSpeed = NULL;
  DO_Set* setEpot = NULL;
  FILE* fp;
  vvTokenizer::TokenType ttype;
  const char* path;
  char* filename;
  char buf[1024];
  int retVal;
  int i, c;
  int iNumber, iType, iMass, iX, iY, iZ, iVX, iVY, iVZ, iEpot;      // locations of numbers in row
  float fNumber, fType, fMass, fX, fY, fZ, fVX, fVY, fVZ, fEpot;    // temporary values 
  int timesteps = 0;
  int atoms;
  float speed;
  bool constrainSpeed;
  float speedLimit[2];  // min and max
  float speedFound[2];  // min and max
  bool periodic;
  float splitter[3];    // splitter for periodic boundaries
  float boundaries[3] = {-1.0f, -1.0f, -1.0f};  // size of box in which atoms are located
  float typeOffset;
  float value;
  int discarded;    // discarded atoms due to speed constraints
  bool ignoreWarnings;

  // Initialize parameters:
  constrainSpeed = pboConstrainSpeed->getValue();
  speedLimit[0] = pfsSpeedMin->getValue();
  speedLimit[1] = pfsSpeedMax->getValue();
  periodic = pboPeriodic->getValue();
  for (i=0; i<3; ++i)
    splitter[i] = pfvPeriodic->getValue(i);
  typeOffset = pfsTypeOffset->getValue();
  ignoreWarnings = pboIgnoreWarnings->getValue();
  speedFound[0] = FLT_MAX;
  speedFound[1] = -FLT_MAX;

  // Open first checkpoint file:
  path = pbrCheckpointFile->getValue();

  if (!vvToolshed::isFile(path))
  {
    sprintf(buf, "Checkpoint file %s not found.", path);
    send_error(buf);
    return STOP_PIPELINE;
  }    

  // Create temporary filename that can be modified to increase:
  filename = new char[strlen(path) + 1];
  strcpy(filename, path);
  
  // Read time steps one by one:
  while (fp=fopen(filename, "rb"))
  {
    // Initialize tokenizer:
    vvTokenizer* tokenizer = new vvTokenizer(fp);
    tokenizer->setEOLisSignificant(true);
    tokenizer->setCaseConversion(vvTokenizer::VV_UPPER);
    tokenizer->setParseNumbers(true);

    // Parse header:
    i = 0;
    iNumber = iType = iMass = iX = iY = iZ = iVX = iVY = iVZ = iEpot = -1;
    while ((ttype = tokenizer->nextToken()) != vvTokenizer::VV_EOF)
    {
      if (strcmp(tokenizer->sval, "#")==0)
      {
        if (!ignoreWarnings)
        {
          sprintf(buf, "Warning: deprecated header type in line %d of checkpoint file %s.", tokenizer->getLineNumber(), filename);
          send_info(buf);
        }
      }
      else if ((strcmp(tokenizer->sval, "#C")==0 || strcmp(tokenizer->sval, "CONTENTS")==0) && i==0)
      {
        // Read data format description:
        while ((ttype = tokenizer->nextToken()) != vvTokenizer::VV_EOL && ttype != vvTokenizer::VV_EOF)
        {
          if (strcmp(tokenizer->sval, "NUMBER")==0)    iNumber = i;
          else if (strcmp(tokenizer->sval, "TYPE")==0) iType = i;
          else if (strcmp(tokenizer->sval, "MASS")==0) iMass = i;
          else if (strcmp(tokenizer->sval, "X")==0)    iX = i;
          else if (strcmp(tokenizer->sval, "Y")==0)    iY = i;
          else if (strcmp(tokenizer->sval, "Z")==0)    iZ = i;
          else if (strcmp(tokenizer->sval, "VX")==0)   iVX = i;
          else if (strcmp(tokenizer->sval, "VY")==0)   iVY = i;
          else if (strcmp(tokenizer->sval, "VZ")==0)   iVZ = i;
          else if (strcmp(tokenizer->sval, "EPOT")==0) iEpot = i;
          ++i;
        }
      }
      else if (strcmp(tokenizer->sval, "#X")==0 || strcmp(tokenizer->sval, "BOX_X")==0 )
      {
        if ((ttype = tokenizer->nextToken()) != vvTokenizer::VV_EOL && ttype != vvTokenizer::VV_EOF)
        {
          boundaries[0] = tokenizer->nval;
        }        
      }
      else if (strcmp(tokenizer->sval, "#Y")==0 || strcmp(tokenizer->sval, "BOX_Y")==0 )
      {
        for (i=0; i<2 && ((ttype = tokenizer->nextToken()) != vvTokenizer::VV_EOL && ttype != vvTokenizer::VV_EOF); ++i)
        {
          boundaries[1] = tokenizer->nval;
        }        
      }
      else if (strcmp(tokenizer->sval, "#Z")==0 || strcmp(tokenizer->sval, "BOX_Z")==0 )
      {
        for (i=0; i<3 && ((ttype = tokenizer->nextToken()) != vvTokenizer::VV_EOL && ttype != vvTokenizer::VV_EOF); ++i)
        {
          boundaries[2] = tokenizer->nval;
        }        
      }
      else if (strcmp(tokenizer->sval, "#E")==0 || strcmp(tokenizer->sval, "ENDHEADER")==0) 
      {
        tokenizer->nextLine();
        break;
      }
      else
      {
        tokenizer->nextLine();      // ignore these types of information lines
      }
    }

    if (iNumber<0 && iType<0 && iMass<0 && iX<0 && iY<0 && iZ<0 && 
        iVX<0 && iVY<0 && iVZ<0 && iEpot<0)
    {
      sprintf(buf, "Data description header line missing in file %s.", filename);
      send_error(buf);
      break;
    }

    if (boundaries[0]<0.0f && boundaries[1]<0.0f && boundaries[2]<0.0f) 
    {
      sprintf(buf, "Boundary box information missing in file %s.", filename);
      send_error(buf);
      break;
    }

    // Parse data area:
    c = 0;
    discarded = 0;
    for(;;)
    {
      ttype = tokenizer->nextToken();
      if (ttype == vvTokenizer::VV_EOL || ttype == vvTokenizer::VV_EOF)
      {
        if (c>0)
        {
          if (constrainSpeed) 
          {
            speed = absVector(fVX, fVY, fVZ);
            if (speed < speedFound[0]) speedFound[0] = speed;
            if (speed > speedFound[1]) speedFound[1] = speed;
          }
          if (!constrainSpeed || (speed >= speedLimit[0] && speed <= speedLimit[1]))   // is atom in valid range?
          {
            number.append(fNumber);
            type.append(fType + typeOffset);
            mass.append(fMass);
            x.append(fX);
            y.append(fY);
            z.append(fZ);
            vx.append(fVX);
            vy.append(fVY);
            vz.append(fVZ);
            epot.append(fEpot);
          }
          else ++discarded;
        }
        if (ttype == vvTokenizer::VV_EOF) break;
        c = 0;
        continue;
      }
      if (ttype != vvTokenizer::VV_NUMBER) 
      {
        sprintf(buf, "Error: cannot parse line %d of file %s.", tokenizer->getLineNumber(), filename);
        send_info(buf);
        break;
      }
      if (c==0)
      {
        fNumber = fType = fMass = fX = fY = fZ = fVX = fVY = fVZ = fEpot = 0.0f;
      }
      if (c==iNumber) fNumber = tokenizer->nval;
      else if (c==iType) fType = tokenizer->nval;
      else if (c==iMass) fMass = tokenizer->nval;
      else if (c==iEpot) fEpot = tokenizer->nval;
      else if (c==iX) 
      {
        if (tokenizer->nval < 0.0f || tokenizer->nval > boundaries[0])
        {
          if (!ignoreWarnings)
          {
            sprintf(buf, "Warning: x coordinate %f out of range in line %d of file %s.", tokenizer->nval, tokenizer->getLineNumber(), filename);
            send_info(buf);
          }
        }
        value = tokenizer->nval;
        if (periodic && value<splitter[0]) value += boundaries[0];
        fX = value;
      }
      else if (c==iY) 
      {
        if (tokenizer->nval < 0.0f || tokenizer->nval > boundaries[1])
        {
          if (!ignoreWarnings)
          {
            sprintf(buf, "Warning: y coordinate %f out of range in line %d of file %s.", tokenizer->nval, tokenizer->getLineNumber(), filename);
            send_info(buf);
          }
        }
        value = tokenizer->nval;
        if (periodic && value<splitter[1]) value += boundaries[1];
        fY = value;
      }
      else if (c==iZ) 
      {
        if (tokenizer->nval < 0.0f || tokenizer->nval > boundaries[2])
        {
          if (!ignoreWarnings)
          {
            sprintf(buf, "Warning: z coordinate %f out of range in line %d of file %s.", tokenizer->nval, tokenizer->getLineNumber(), filename);
            send_info(buf);
          }
        }
        value = tokenizer->nval;
        if (periodic && value<splitter[2]) value += boundaries[2];
        fZ = value;
      }
      else if (c==iVX) fVX = tokenizer->nval;
      else if (c==iVY) fVY = tokenizer->nval;
      else if (c==iVZ) fVZ = tokenizer->nval;
      ++c;  
    }
    delete tokenizer;    
    fclose(fp);
    
    // Create Covise data objects from arrays:
    atoms = number.count();
    ++timesteps;
    
    sprintf(buf, "%s_%d", poType->getObjName(), timesteps);
    doType = new DO_Unstructured_S3D_Data(buf, atoms, type.getArrayPtr());
    aType.append(doType);

    sprintf(buf, "%s_%d", poMass->getObjName(), timesteps);
    doMass = new DO_Unstructured_S3D_Data(buf, atoms, mass.getArrayPtr());
    aMass.append(doMass);

    sprintf(buf, "%s_%d", poPoints->getObjName(), timesteps);
    doPoints = new DO_Points(buf, atoms, x.getArrayPtr(), y.getArrayPtr(), z.getArrayPtr());
    aPoints.append(doPoints);

    sprintf(buf, "%s_%d", poSpeed->getObjName(), timesteps);
    doSpeed = new DO_Unstructured_V3D_Data(buf, atoms, vx.getArrayPtr(), vy.getArrayPtr(), vz.getArrayPtr());
    aSpeed.append(doSpeed);

    sprintf(buf, "%s_%d", poEpot->getObjName(), timesteps);
    doEpot = new DO_Unstructured_S3D_Data(buf, atoms, epot.getArrayPtr());
    aEpot.append(doEpot);

    // Clear raw data arrays:
    number.clear();
    type.clear();
    mass.clear();
    x.clear();
    y.clear();
    z.clear();
    vx.clear();
    vy.clear();
    vz.clear();
    epot.clear();

    // Print info message:
    sprintf(buf, "%d atoms loaded from checkpoint file %s.", atoms, filename);
    send_info(buf);
    if (constrainSpeed)
    {
      sprintf(buf, "%d atoms discarded due to speed constraints.", discarded);
      send_info(buf);
    }
    sprintf(buf, "Size of simulation box: %f x %f x %f", boundaries[0], boundaries[1], boundaries[2]);
    send_info(buf);

    // Process next time step:
    if (!vvToolshed::increaseFilename(filename)) break;
  }

  // Terminate data object arrays:
  aType.append(NULL);
  aMass.append(NULL);
  aPoints.append(NULL);
  aSpeed.append(NULL);
  aEpot.append(NULL);

  if (constrainSpeed)
  {
    sprintf(buf, "Minimum speed: %f, Maximum speed: %f", speedFound[0], speedFound[1]);
    send_info(buf);
  }

  if (timesteps==0)
  {
    send_error("No atoms loaded.");
    retVal = STOP_PIPELINE;
  }    
  else    // data has been loaded and can now be converted to sets
  {
    // Create set objects:
    setPoints = new DO_Set(poPoints->getObjName(), (DistributedObject**)aPoints.getArrayPtr());
    setType   = new DO_Set(poType->getObjName(),   (DistributedObject**)aType.getArrayPtr());
    setMass   = new DO_Set(poMass->getObjName(),   (DistributedObject**)aMass.getArrayPtr());
    setSpeed  = new DO_Set(poSpeed->getObjName(),  (DistributedObject**)aSpeed.getArrayPtr());
    setEpot   = new DO_Set(poEpot->getObjName(),   (DistributedObject**)aEpot.getArrayPtr());

    // Now the arrays can be cleared:
    aPoints.clear();
    aType.clear();
    aMass.clear();
    aSpeed.clear();
    aEpot.clear();

    // Set timestep attribute:
    if (timesteps > 1) 
    {
        sprintf(buf, "%d %d", 0, timesteps-1);
        setPoints->set_attribute("TIMESTEP", buf);
        setType->set_attribute("TIMESTEP", buf);
        setMass->set_attribute("TIMESTEP", buf);
        setSpeed->set_attribute("TIMESTEP", buf);
        setEpot->set_attribute("TIMESTEP", buf);
    }

    // Assign sets to output ports:
    poPoints->setObj(setPoints);
    poType->setObj(setType);
    poMass->setObj(setMass);
    poSpeed->setObj(setSpeed);
    poEpot->setObj(setEpot);

    sprintf(buf, "Timesteps loaded: %d", timesteps);
    send_info(buf);

    retVal = CONTINUE_PIPELINE;
  }
  delete[] filename;
  return retVal;
}

