/**************************************************************************\
**                                                           (C)2002 RUS  **
**                                                                        **
** Description: Read IMD checkpoint files from ITAP.                      **
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
coCheckpointFile::coCheckpointFile(coModule* m, FILE* file)
{ 
  int i,j;
  for (i=0; i<3; ++i) 
    for (j=0; j<3; ++j)
      boxSize[i][j] = -1.0f;
  module = m;
  fp = file;
}

/// Destructor
coCheckpointFile::~coCheckpointFile()
{
  clearParamNames();
}

/// Remove parameter names array from memory
void coCheckpointFile::clearParamNames()
{
  for (int i=0; i<paramNames.count(); ++i)
    delete[] paramNames[i];
  paramNames.clear();
}

/** @return index of passed parameter string (>=0) or -1 if text was not found
*/
int coCheckpointFile::getParamIndex(char* text)
{
  for (int i=0; i<paramNames.count(); ++i)
  {
    if (strcmp(paramNames[i], text) == 0)
      return i;
  }
  return -1;
}

/// @return index of speed components or -1 if not in paramNames
void coCheckpointFile::getSpeedIndices(int* vx, int* vy, int* vz)
{
  *vx = *vy = *vz = -1;

  for (int i=0; i<paramNames.count(); ++i)
  {
    if (strcmp(paramNames[i], "vx") == 0)
      *vx = i;
    else if (strcmp(paramNames[i], "vy") == 0)
      *vy = i;
    else if (strcmp(paramNames[i], "vz") == 0)
      *vz = i;
  }
}

/// @return index of location components or -1 if not in paramNames
void coCheckpointFile::getLocationIndices(int* x, int* y, int* z)
{
  *x = *y = *z = -1;

  for (int i=0; i<paramNames.count(); ++i)
  {
    if (strcmp(paramNames[i], "x") == 0)
      *x = i;
    else if (strcmp(paramNames[i], "y") == 0)
      *y = i;
    else if (strcmp(paramNames[i], "z") == 0)
      *z = i;
  }
}

/** Parse file header and move file pointer to beginning of data area.
  @return true if successful, false on error
*/
bool coCheckpointFile::parseHeader()
{
  vvTokenizer::TokenType ttype;
  char buf[1024];
  char* newParam;
  int i,j;

  if (fp==NULL) return false;

  clearParamNames();
  fseek(fp, 0, SEEK_SET);

  // Initialize tokenizer:
  vvTokenizer* tokenizer = new vvTokenizer(fp);
  tokenizer->setEOLisSignificant(true);
  tokenizer->setCaseConversion(vvTokenizer::VV_LOWER);
  tokenizer->setParseNumbers(true);

  // Initialize box size:
  for (i=0; i<3; ++i) 
    for (j=0; j<3; ++j)
      boxSize[i][j] = -1.0f;

  // Parse header:
  while ((ttype = tokenizer->nextToken()) != vvTokenizer::VV_EOF)
  {
    if (strcmp(tokenizer->sval, "#")==0)
    {
      sprintf(buf, "Warning: invalid checkpoint header type in line %d.", tokenizer->getLineNumber());
      module->send_info(buf);
    }
    else if ((strcmp(tokenizer->sval, "#c")==0 || strcmp(tokenizer->sval, "contents")==0))
    {
      // Read data format description:
      while ((ttype = tokenizer->nextToken()) != vvTokenizer::VV_EOL && ttype != vvTokenizer::VV_EOF)
      {
        newParam = new char[strlen(tokenizer->sval) + 1];
        strcpy(newParam, tokenizer->sval);
        paramNames.append(newParam);
      }
    }
    else if (strcmp(tokenizer->sval, "#x")==0 || strcmp(tokenizer->sval, "box_x")==0 )
    {
      for (i=0; i<3; ++i)
      {
        ttype = tokenizer->nextToken();
        if (ttype == vvTokenizer::VV_EOL || ttype == vvTokenizer::VV_EOF) break;
        boxSize[0][i] = tokenizer->nval;
      }        
    }
    else if (strcmp(tokenizer->sval, "#y")==0 || strcmp(tokenizer->sval, "box_y")==0 )
    {
      for (i=0; i<3; ++i)
      {
        ttype = tokenizer->nextToken();
        if (ttype == vvTokenizer::VV_EOL || ttype == vvTokenizer::VV_EOF) break;
        boxSize[1][i] = tokenizer->nval;
      }        
    }
    else if (strcmp(tokenizer->sval, "#z")==0 || strcmp(tokenizer->sval, "box_z")==0 )
    {
      for (i=0; i<3; ++i)
      {
        ttype = tokenizer->nextToken();
        if (ttype == vvTokenizer::VV_EOL || ttype == vvTokenizer::VV_EOF) break;
        boxSize[2][i] = tokenizer->nval;
      }        
    }
    else if (strcmp(tokenizer->sval, "#e")==0 || strcmp(tokenizer->sval, "endheader")==0) 
    {
      tokenizer->nextLine();
      break;
    }
    else if (ttype != vvTokenizer::VV_EOL)    // don't ignore next line if we are already at the end of the line
    {
      tokenizer->nextLine();      // ignore the rest of the line
    }
  }

  delete tokenizer;
  tokenizer = NULL;

  // Verify bounding box elements:
  for (i=0; i<3; ++i)
    for (j=0; j<3; ++j)
    {
      if (boxSize[i][j] != 0.0f && i!=j)
      {
        if (((coReadIMD*)module)->displayWarnings())
        {
          sprintf(buf, "Ignoring non-zero bounding box parameter not on diagonal: %f", boxSize[i][j]);
          module->send_info(buf);
        }
      }
      else if (i==j && boxSize[i][j]<=0.0f)
      {
        if (((coReadIMD*)module)->displayWarnings())
        {
          sprintf(buf, "Non-positive bounding box parameter on diagonal set to 1.0: %f", boxSize[i][j]);
          module->send_info(buf);
        }
        boxSize[i][j] = 1.0f;
      }
    }
    
  return true;
}

/// Constructor
coReadIMD::coReadIMD()
{
  char* defParam[] = {"none"};

  set_module_description("Read IMD checkpoint files to create a list of points and scalar parameters.");

  // Create ports:
  poPoints = addOutPort("Location", "Set_Points", "Atom location");
  poPoints->setInfo("Atom location");

  poSpeed = addOutPort("Speed", "Set_Unstructured_V3D_Data", "Atom speed");
  poSpeed->setInfo("Atom speed");

  poScalar1 = addOutPort("Scalar1", "Set_Unstructured_S3D_Data", "Scalar parameter #1");
  poScalar1->setInfo("Scalar parameter #1");

  poScalar2 = addOutPort("Scalar2", "Set_Unstructured_S3D_Data", "Scalar parameter #2");
  poScalar2->setInfo("Scalar parameter #2");

  // Create parameters:
  pbrCheckpointFile = addBrowserParam("FilePath", "First IMD checkpoint file of sequence or single file");
  pbrCheckpointFile->setValue("data/","*.cpt");
  pbrCheckpointFile->setImmediate(true);

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

  pchScalar1 = addChoiceParam("Parameter1", "Parameter at first scalar output port");
  pchScalar1->setValue(1, defParam, 1);

  pchScalar2 = addChoiceParam("Parameter2", "Parameter at second scalar output port");
  pchScalar2->setValue(1, defParam, 1);
  
  pfsFactor1 = addFloatScalarParam("Factor1", "Factor to multiply scalar parameter #1 with");
  pfsFactor1->setValue(1.0f);

  pfsFactor2 = addFloatScalarParam("Factor2", "Factor to multiply scalar parameter #2 with");
  pfsFactor2->setValue(1.0f);

  pfsOffset1 = addFloatScalarParam("Offset1", "Offset to add to scalar parameter #1");
  pfsOffset1->setValue(0.0f);

  pfsOffset2 = addFloatScalarParam("Offset2", "Offset to add to scalar parameter #2");
  pfsOffset2->setValue(0.0f);

  pboWarnings = addBooleanParam("Warnings", "Display warnings when reading files");
  pboWarnings->setValue(true);
}

void coReadIMD::param(const char* paramName)
{
  if (strcmp(paramName, pbrCheckpointFile->getName())==0)
  {
    char* defaultParam[2];
    defaultParam[0] = new char[strlen(paramNames[pchScalar1->getValue() - 1]) + 1];
    defaultParam[1] = new char[strlen(paramNames[pchScalar2->getValue() - 1]) + 1];
    strcpy(defaultParam[0], paramNames[pchScalar1->getValue() - 1]);
    strcpy(defaultParam[1], paramNames[pchScalar2->getValue() - 1]);
    readParameters();                 // parse header for names of atom parameters
    updateParameters(defaultParam);   // update parameters in module's choice lists
    delete[] defaultParam[1];
    delete[] defaultParam[0];
  }
}

/// Find out which parameters the selected file contains.
void coReadIMD::readParameters()
{
  FILE* fp;
  coCheckpointFile* cpFile;
  char* newString;
  char buf[1024];
  int i;

  fp = fopen(pbrCheckpointFile->getValue(), "rb");
  if (fp==NULL)
  {
    sprintf(buf, "Cannot open source file: %s", pbrCheckpointFile->getValue());
    send_error(buf);
    return;
  }  
  cpFile = new coCheckpointFile(this, fp);
  cpFile->parseHeader();

  // Copy parameter names from checkpoint file class to this class:
  paramNames.clear();
  for (i=0; i<cpFile->paramNames.count(); ++i)
  {
    newString = new char[strlen(cpFile->paramNames[i]) + 1];
    strcpy(newString, cpFile->paramNames[i]);
    paramNames.append(newString);
  }
  delete cpFile;
  cpFile = NULL;
  fclose(fp);
}

/** Update parameter lists in choice boxes.
  @param defaultParam default parameter names: use if they occur in list, 
                      otherwise default to first parameter in list
*/
void coReadIMD::updateParameters(char** defaultParam)
{
  bool found[2] = {false, false};
  int  defIndex[2] = {1, 1};
  for (int i=0; i<paramNames.count(); ++i)
  {
    if (!found[0] && strcmp(paramNames[i], defaultParam[0])==0) 
    {
      defIndex[0] = i + 1;
      found[0] = true;
    }
    if (!found[1] && strcmp(paramNames[i], defaultParam[1])==0) 
    {
      defIndex[1] = i + 1;
      found[1] = true;
    }
  }
  pchScalar1->setValue(paramNames.count(), paramNames.getArrayPtr(), defIndex[0]);
  pchScalar2->setValue(paramNames.count(), paramNames.getArrayPtr(), defIndex[1]);
}

/// @return absolute value of a vector
float coReadIMD::absVector(float x, float y, float z)
{
  return sqrt(x * x + y * y + z * z);
}

/// @return true if warnings are to be displayed
bool coReadIMD::displayWarnings()
{
  return pboWarnings->getValue();
}

/// Compute routine: load checkpoint file
int coReadIMD::compute()
{
  const int DEF_DIM = 10000;
  vvArray<DO_Points*> aPoints(1,5);
  vvArray<DO_Unstructured_V3D_Data*> aSpeed(1,5);
  vvArray<DO_Unstructured_S3D_Data*> aScalar1(1,5);
  vvArray<DO_Unstructured_S3D_Data*> aScalar2(1,5);
  vvArray<float> x(0,DEF_DIM);
  vvArray<float> y(0,DEF_DIM);
  vvArray<float> z(0,DEF_DIM);
  vvArray<float> vx(0,DEF_DIM);
  vvArray<float> vy(0,DEF_DIM);
  vvArray<float> vz(0,DEF_DIM);
  vvArray<float> scalar1(0,DEF_DIM);
  vvArray<float> scalar2(0,DEF_DIM);
  DO_Points* doPoints = NULL;
  DO_Unstructured_V3D_Data* doSpeed = NULL;
  DO_Unstructured_S3D_Data* doScalar1 = NULL;
  DO_Unstructured_S3D_Data* doScalar2 = NULL;
  DO_Set* setPoints = NULL;
  DO_Set* setSpeed = NULL;
  DO_Set* setScalar1 = NULL;
  DO_Set* setScalar2 = NULL;
  FILE* fp;
  const char* path;
  char* filename;
  char buf[1024];
  int retVal;
  int i, c;
  int iX, iY, iZ, iVX, iVY, iVZ, iScalar[2];      // locations of numbers in row
  float fX, fY, fZ, fVX, fVY, fVZ, fScalar[2];    // temporary values 
  int timesteps = 0;
  int atoms;
  float speed;
  bool constrainSpeed;
  float speedLimit[2];  // min and max
  float speedFound[2];  // min and max
  bool periodic;
  float splitter[3];    // splitter for periodic boundaries
  float value;
  int discarded;    // discarded atoms due to speed constraints
  bool warnings;
  coCheckpointFile* cpFile;
  vvTokenizer::TokenType ttype;

  // Initialize parameters:
  constrainSpeed = pboConstrainSpeed->getValue();
  speedLimit[0] = pfsSpeedMin->getValue();
  speedLimit[1] = pfsSpeedMax->getValue();
  periodic = pboPeriodic->getValue();
  for (i=0; i<3; ++i)
    splitter[i] = pfvPeriodic->getValue(i);
  warnings = pboWarnings->getValue();
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
    cpFile = new coCheckpointFile(this, fp);
    cpFile->parseHeader();
    iScalar[0] = cpFile->getParamIndex(paramNames[pchScalar1->getValue()-1]);
    iScalar[1] = cpFile->getParamIndex(paramNames[pchScalar2->getValue()-1]);
    cpFile->getLocationIndices(&iX, &iY, &iZ);
    if (iX<0 || iY<0 || iZ<0) 
    {
      send_info("Atom position missing in header.");
      delete cpFile;
      break;
    }
    cpFile->getSpeedIndices(&iVX, &iVY, &iVZ);

    // Initialize tokenizer:
    vvTokenizer* tokenizer = new vvTokenizer(fp);
    tokenizer->setEOLisSignificant(true);
    tokenizer->setCaseConversion(vvTokenizer::VV_LOWER);
    tokenizer->setParseNumbers(true);

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
            x.append(fX);
            y.append(fY);
            z.append(fZ);
            vx.append(fVX);
            vy.append(fVY);
            vz.append(fVZ);
            if (iScalar[0]>-1) scalar1.append(fScalar[0] * pfsFactor1->getValue() + pfsOffset1->getValue());
            if (iScalar[1]>-1) scalar2.append(fScalar[1] * pfsFactor2->getValue() + pfsOffset2->getValue());
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
        fX = fY = fZ = fVX = fVY = fVZ = fScalar[0] = fScalar[1] = 0.0f;
      }

      // Check for user defined parameters:
      if (c==iScalar[0]) fScalar[0] = tokenizer->nval;
      if (c==iScalar[1]) fScalar[1] = tokenizer->nval;

      // Check for position and velocity:
      if (c==iX) 
      {
        if (tokenizer->nval < 0.0f || tokenizer->nval > cpFile->boxSize[0][0])
        {
          if (warnings)
          {
            sprintf(buf, "Warning: x coordinate %f out of range in line %d of file %s.", 
              tokenizer->nval, tokenizer->getLineNumber(), filename);
            send_info(buf);
          }
        }
        value = tokenizer->nval;
        if (periodic && value<splitter[0]) value += cpFile->boxSize[0][0];
        fX = value;
      }
      else if (c==iY) 
      {
        if (tokenizer->nval < 0.0f || tokenizer->nval > cpFile->boxSize[1][1])
        {
          if (warnings)
          {
            sprintf(buf, "Warning: y coordinate %f out of range in line %d of file %s.", 
              tokenizer->nval, tokenizer->getLineNumber(), filename);
            send_info(buf);
          }
        }
        value = tokenizer->nval;
        if (periodic && value<splitter[1]) value += cpFile->boxSize[1][1];
        fY = value;
      }
      else if (c==iZ) 
      {
        if (tokenizer->nval < 0.0f || tokenizer->nval > cpFile->boxSize[2][2])
        {
          if (warnings)
          {
            sprintf(buf, "Warning: z coordinate %f out of range in line %d of file %s.", tokenizer->nval, tokenizer->getLineNumber(), filename);
            send_info(buf);
          }
        }
        value = tokenizer->nval;
        if (periodic && value<splitter[2]) value += cpFile->boxSize[2][2];
        fZ = value;
      }
      else if (c==iVX) fVX = tokenizer->nval;
      else if (c==iVY) fVY = tokenizer->nval;
      else if (c==iVZ) fVZ = tokenizer->nval;
      ++c;  
    }
    delete tokenizer;    
    tokenizer = NULL;
    fclose(fp);
    
    // Create Covise data objects from arrays:
    atoms = x.count();
    if (atoms != y.count() || atoms != z.count()) 
    { 
      send_error("Error in coordinate list.");
      break;
    }
    ++timesteps;
    
    sprintf(buf, "%s_%d", poPoints->getObjName(), timesteps);
    doPoints = new DO_Points(buf, atoms, x.getArrayPtr(), y.getArrayPtr(), z.getArrayPtr());
    aPoints.append(doPoints);

    sprintf(buf, "%s_%d", poSpeed->getObjName(), timesteps);
    doSpeed = new DO_Unstructured_V3D_Data(buf, atoms, vx.getArrayPtr(), vy.getArrayPtr(), vz.getArrayPtr());
    aSpeed.append(doSpeed);

    sprintf(buf, "%s_%d", poScalar1->getObjName(), timesteps);
    doScalar1 = new DO_Unstructured_S3D_Data(buf, atoms, scalar1.getArrayPtr());
    aScalar1.append(doScalar1);

    sprintf(buf, "%s_%d", poScalar2->getObjName(), timesteps);
    doScalar2 = new DO_Unstructured_S3D_Data(buf, atoms, scalar2.getArrayPtr());
    aScalar2.append(doScalar2);

    // Clear raw data arrays:
    x.clear();
    y.clear();
    z.clear();
    vx.clear();
    vy.clear();
    vz.clear();
    scalar1.clear();
    scalar2.clear();

    // Print info message:
    sprintf(buf, "%d atoms loaded from checkpoint file %s.", atoms, filename);
    send_info(buf);
    if (constrainSpeed)
    {
      sprintf(buf, "%d atoms discarded due to speed constraints.", discarded);
      send_info(buf);
    }
    sprintf(buf, "Size of simulation box: %f x %f x %f", 
      cpFile->boxSize[0][0], cpFile->boxSize[1][1], cpFile->boxSize[2][2]);
    send_info(buf);

    delete cpFile; cpFile = NULL;

    // Process next time step:
    if (!vvToolshed::increaseFilename(filename)) break;
  }

  // Terminate data object arrays:
  aPoints.append(NULL);
  aSpeed.append(NULL);
  aScalar1.append(NULL);
  aScalar2.append(NULL);

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
    setPoints  = new DO_Set(poPoints->getObjName(),  (DistributedObject**)aPoints.getArrayPtr());
    setSpeed   = new DO_Set(poSpeed->getObjName(),   (DistributedObject**)aSpeed.getArrayPtr());
    setScalar1 = new DO_Set(poScalar1->getObjName(), (DistributedObject**)aScalar1.getArrayPtr());
    setScalar2 = new DO_Set(poScalar2->getObjName(), (DistributedObject**)aScalar2.getArrayPtr());

    // Now the arrays can be cleared:
    aPoints.clear();
    aSpeed.clear();
    aScalar1.clear();
    aScalar2.clear();

    // Set timestep attribute:
    if (timesteps > 1) 
    {
        sprintf(buf, "%d %d", 0, timesteps-1);
        setPoints->set_attribute("TIMESTEP", buf);
        setSpeed->set_attribute("TIMESTEP", buf);
        setScalar1->set_attribute("TIMESTEP", buf);
        setScalar2->set_attribute("TIMESTEP", buf);
    }

    // Assign sets to output ports:
    poPoints->setObj(setPoints);
    poSpeed->setObj(setSpeed);
    poScalar1->setObj(setScalar1);
    poScalar2->setObj(setScalar2);

    sprintf(buf, "Timesteps loaded: %d", timesteps);
    send_info(buf);

    retVal = CONTINUE_PIPELINE;
  }
  delete[] filename;
  return retVal;
}

