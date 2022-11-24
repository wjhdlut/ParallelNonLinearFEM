/* ----------------------------------------------------------------------
 *| This file createS a class to read the input files and output the     |
 *| results into files. Then the  of displacement, stress and strain     |
 *| can be pbtained by software Tecplot.                                 |
 *|                                                                      |
 *| - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
 *| Filename: BoundaryCondition.cpp                                      |
 *|   Author: WangJianHua                                                |
 *|     Date: 2022-3-20.                                                 |
 *| All rights reserved!                                                 |
  ----------------------------------------------------------------------- */
#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

#include<mpi.h>
#include<string>
#include<vector>
#include<unordered_map>

#include"ElementLib.h"

using namespace std;

#define MAX_LENGTH_CHAR 1000

struct node{
    int ID;
    vector<double> coord;
};

struct nodeBoundaryCondition{
    string Dof;
    vector<double> value;
};

struct elemDistributeForce
{
    int nodeID;
    vector<double> value;
};

class InputOutput
{
public:
    explicit InputOutput(int rank);
    ~InputOutput();

public:
    int GetNodeNum();

    int GetElemNum();

    int GetAnalyseType();

    int GetAxisOfSymmetry();

    int GetLargeStrain();

    int GetSolAlgorithmType();

    int GetSolArcLengthType();

    unordered_map<int, vector<double>> GetNodeCoord();

    unordered_map<int, vector<int>> GetElemType();

    unordered_map<int, vector<int>> GetElemGroup();

    unordered_map<int, nodeBoundaryCondition> GetPrescribeNodeDisp();

    unordered_map<int, nodeBoundaryCondition> GetPrescribeNodeForce();

    unordered_map<int, vector<int>> GetElemConnect();

    unordered_map<int, vector<elemDistributeForce>> GetDistributeForce();

    vector<vector<double>> GetIncrements();

public:
    void ReadInput();

private:
    string  m_inputFileName;
    string  *m_outputFileName;
    int     m_rank;
    int     m_argv;
    char    **m_args;

    int     m_analyseType;
    int     m_axisOfSymmetry;
    int     m_largeStrain;
    int     m_solAlgorithmType;
    int     m_solArcLengthType;
    int     m_numElemGroup;
    int     m_numElemType;
    int     m_numElem;
    int     m_numNode;
    int     m_numPrescribedDispNode;
    int     m_numPrescribedForceNode;
    int     m_numMat;
    int     m_numLoadStep;

    bool    noTitle;
    bool    noSolAlgorithmType;
    bool    noSolArcLengthType;
    bool    cylindrical;
    bool    cartes;

private:
    unordered_map<int, vector<int>>                 m_elemType;
    unordered_map<int, vector<double>>              m_nodeSet;
    unordered_map<int, vector<int>>                 m_elemGroup;
    unordered_map<int, nodeBoundaryCondition>       m_prescribeNodeDisp;
    unordered_map<int, nodeBoundaryCondition>       m_prescribeNodeForce;
    unordered_map<int, vector<int>>                 m_elemConnect;
    unordered_map<int, vector<elemDistributeForce>> m_distributeForce;
    vector<vector<double>>                          m_increments;

private:
    /**
     * @brief PrintBasicInfo
     */
    void PrintBasicInfo();
    /**
     * @brief Set Input File Name
     */
    void SetInputFileName();

    /**
     * @brief InitialVariables
     */
    void InitialVariables();

    /**
     * @brief Read Input file information
     * @param lines
     */
    void ReadInputString(vector<string>&lines);

    /**
     * @brief ReadTitleInfo
     * @param strData
     */
    void ReadTitleInfo(const string&strData);

    /**
     * @brief Read Analyse Type Information
     * @param strData
     */
    void ReadAnalyseTypeInfo(const string&strData);

    /**
     * @brief Read Axis Of Symmetry
     * @param strData
     */
    void ReadAxisOfSymmetryInfo(const string&strData);

    /**
     * @brief Read Large Strain Information
     * @param strData
     */
    void ReadLargeStrainInfo(const string&strData);

    /**
     * @brief ReadSolAlgorithmInfo
     * @param strData
     */
    void ReadSolAlgorithmInfo(const string&strData);

    /**
     * @brief ReadArcLengthInfo
     * @param strData
     */
    void ReadArcLengthInfo(const string&strData);

    /**
     * @brief ReadElemGroupInfo
     * @param strData
     */
    void ReadElemGroupInfo(const string&strData);

    /**
     * @brief ReadElemTypeInfo
     * @param strData
     */
    void ReadElemTypeInfo(const string&strData);

    /**
     * @brief ReadElemConnectInfo
     * @param strData
     */
    void ReadElemConnectInfo(const string&strData);

    /**
     * @brief ReadNodeCoord
     * @param strData
     */
    void ReadNodeCoord(ifstream &nodeFile);

    /**
     * @brief ReadDispBoundaryInfo
     * @param strData
     */
    void ReadDispBoundaryInfo(const string&strData);

    /**
     * @brief ReadMaterialPara
     * @param strData
     */
    void ReadMaterialInfo(const string&strData);

    /**
     * @brief CheckInputFile
     */
    void CheckInputFile();

    /**
     * @brief RegexPattern
     * @param temp_map
     * @param strData
     */
    void RegexPattern(unordered_map<int, vector<int>>&tempMap, const string &strData);

    void ReadMaterialPara(const string&matType, const string&para);

    void ReadGeometryInfo(const string&inputFile);

    void ReadNewtonRaphsonAlgorithmInfo(const string&strData);
};

#endif // INPUTOUTPUT_H
