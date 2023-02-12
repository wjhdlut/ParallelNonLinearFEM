#ifndef INPUTREADER_H
#define INPUTREADER_H

#include <util/DataStructure.h>

namespace NONLINEARFEMIO
{
/**
 * @Brief: Read input file
 * 
 * @param rank 
 * @return GlobalData* 
 */
GlobalData* InputReader(int rank);
};

#endif // INPUTREADER_H

