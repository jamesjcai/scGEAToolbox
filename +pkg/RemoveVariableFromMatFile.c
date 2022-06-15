#include "mex.h"
#include "mat.h"

/* This function removes one or more variables from a MAT file
 * Compile it with
 * >>mex rmvarMatfileMEX.c
 * Afterwards call it with
 * >> rmvarMatfileMEX(FILENAME_WITH_EXTENSION,...variables....)
 * e.g.
 * >> rmvarMatfileMEX('MyFile.mat','var1','var2',...)
 * https://www.mathworks.com/matlabcentral/answers/254048-how-can-i-delete-variables-in-my-mat-file-without-loading-them-into-matlab-7-2-r2006a
 */ 

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    MATFile *f;
    char *filename;
    char *vname;
    int tmp;

      if (nrhs >= 2 )
      {

              if (!mxIsChar(prhs[0]) || !mxIsChar(prhs[0]))
              {
                  mexErrMsgIdAndTxt("RemoveVariableFromMatFile:ClassInputArguments","This function expects the inputs to be char."); 
              }

              filename = mxArrayToString(prhs[0]);
              f =  matOpen(filename,"u");

              if (f == NULL)
              {
                  mxFree(filename);
                  mexErrMsgIdAndTxt("RemoveVariableFromMatFile:UnableToOpenFile","Could not open file. Make sure the file exists and is accessible.");  
              }

              for (int i=1;i<nrhs;i++)
              {

                      vname = mxArrayToString(prhs[i]);
                      tmp = matDeleteVariable(f,vname);

                      if ( tmp != 0)
                      {
                          mexWarnMsgIdAndTxt("RemoveVariableFromMatFile:UnableToDeleteVariable","Could not delete variable. Make sure that the variable exists.");     
                      }

                      mxFree(vname);
                      vname=NULL;

              }

              matClose(f);
              mxFree(filename);

       }
  }