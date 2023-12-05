/**
 * @File Name:     InputReader.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef INPUTREADER_H
#define INPUTREADER_H

#include "../../util/include/DataStructure.h"

namespace NONLINEARFEMIO
{
/**
 * @Brief: Read input file
 * 
 * @param rank 
 * @return GlobalData* 
 */
GlobalData* InputReader(int rank, char **args);
};

#endif // INPUTREADER_H

