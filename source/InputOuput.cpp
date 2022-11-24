#include<vector>
#include<iostream>
#include<fstream>
#include<regex>
#include<algorithm>
#include<math.h>

#include"../include/InputOuput.h"

InputOutput::InputOutput(int rank):m_rank(rank)
{
    //m_inputFileName = ;
    PrintBasicInfo();
    SetInputFileName();
    InitialVariables();
}

InputOutput::~InputOutput()
{
}

void InputOutput::PrintBasicInfo()
{
    cout << "  ________________________________________________________________________ " << endl;
    cout << " |_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|                                          |" << endl;
    cout << " |_|_|_|_|_|_|_|_|_|_|                                                    |" << endl;
    cout << " |_|_|_|_|_|_|                                                            |" << endl;
    cout << " |_|_|_|_|                                                                |" << endl;
    cout << " |_|_|_|               Parallel Nonlinear FEM Package                     |" << endl;
    cout << " |_|_|                ================================                    |" << endl;
    cout << " |_|_|                                                                   _|" << endl;
    cout << " |_|_|                                                                  |_|" << endl;
    cout << " |_|_|_                                                                 |_|" << endl;
    cout << " |_|_|_|_     SMALL AND LARGE STRAIN FINITE ELEMENT ANALYSI            _|_|" << endl;
    cout << " |_|_|_|_|      OF HYPERELASTIC AND ELASTO-PLASTIC SOLIDS             |_|_|" << endl;
    cout << " |_|_|_|_|_                                                        _ _|_|_|" << endl;
    cout << " |_|_|_|_|_|                                                     _|_|_|_|_|" << endl;
    cout << " |_|_|_|_|_|_ _                                              _ _|_|_|_|_|_|" << endl;
    cout << " |_|_|_|_|_|_|_|____________________________________________|_|_|_|_|_|_|_|" << endl;
    cout << " |                                                                        |" << endl;
    cout << " |    Copyright (c) 1996-2008   EA de Souza Neto, D Peric & DRJ Owen      |" << endl;
    cout << " |________________________________________________________________________|" << endl;
    cout << " |                                                                        |" << endl;
    cout << " |    Companion to the textbook:                                          |" << endl;
    cout << " |    EA de Souza Neto, D Peric & DRJ Owen. Computational Methods for     |" << endl;
    cout << " |    Plasticity: Theory and Applications. Wiley, Chichester, 2008.       |" << endl;
    cout << " |________________________________________________________________________|" << endl;
}

void InputOutput::SetInputFileName()
{
    cout << "\n\n\t" << "   Data file name must have extension .dat or .DAT" << endl;
    cout << "\t\t and must not contain blank spaces." <<endl;
    cout << "\t\t Type EXIT ot QUIT to temianl" << endl;
    cout << "\n\t>>> Please input the filename of your input:>>>";
    //    cin >> m_inputFileName;
    m_inputFileName = "/home/wangjianhua/ExampleCpp/ParallelNonLinearFEM/input/ParallelFEM.in";
    cout << "\n\t The input filename is " << m_inputFileName << endl;
}

void InputOutput::InitialVariables()
{
    m_analyseType            = 0;
    m_axisOfSymmetry         = 0;
    m_largeStrain            = 0;
    m_solAlgorithmType       = 0;
    m_solArcLengthType       = 0;
    m_numElemGroup           = 0;
    m_numElemType            = 0;
    m_numElem                = 0;
    m_numNode                = 0;
    m_numPrescribedDispNode  = 0;
    m_numPrescribedForceNode = 0;
    m_numMat                 = 0;
    m_numLoadStep            = 0;

    noTitle                 = true;
    noSolAlgorithmType      = true;
    noSolArcLengthType      = true;
    cylindrical             = true;
    cartes                  = true;

    m_elemType              = {};
    m_nodeSet               = {};
    m_elemGroup             = {};
    m_prescribeNodeDisp     = {};
    m_prescribeNodeForce    = {};
    m_elemConnect           = {};
    m_distributeForce       = {};
    m_increments            = {};
    m_increments            = {};
}

void InputOutput::ReadInput()
{
    if (m_rank == 0)
    {
        vector<string> lines{};
        ReadInputString(lines);

        for (int var = 0; var < lines.size(); var++) {
            // read title
            if (string::npos != lines[var].find("TITLE")) {
                ReadTitleInfo(lines[var]);
            }

            // read node, mesh, boundary information
            if (string::npos != lines[var].find("INCLUDE")) {
                ReadGeometryInfo(lines[var]);
            }

            // read analysis type
            if (string::npos != lines[var].find("ANALYSIS_TYPE")) {
                ReadAnalyseTypeInfo(lines[var]);
                if (3 == m_analyseType) {
                    ReadAxisOfSymmetryInfo(lines[++var]);
                }
            }

            // read large strain OFF or ON
            if (string::npos != lines[var].find("LARGE_STRAIN_FORMULATION")) {
                ReadLargeStrainInfo(lines[var]);
            }

            // read solution algorithm
            if (string::npos != lines[var].find("SOLUTION_ALGORITHM")) {
                ReadSolAlgorithmInfo(lines[var]);
                if (m_solAlgorithmType < 0) {
                    ReadArcLengthInfo(lines[++var]);
                }
            }

            // element_group
            if (string::npos != lines[var].find("ELEMENT_GROUPS")) {
                ReadElemGroupInfo(lines[var]);
            }

            // read element type
            if (string::npos != lines[var].find("ELEMENT_TYPES")) {
                ReadElemTypeInfo(lines[var]);
            }

            if (string::npos != lines[var].find("MATERIALS")) {
                ReadMaterialInfo(lines[var]);
            }

            if (string::npos != lines[var].find("INCREMENTS")) {
                ReadNewtonRaphsonAlgorithmInfo(lines[var]);
            }

        }

        CheckInputFile();

    }
}

void InputOutput::ReadInputString(vector<string>&lines)
{
    ifstream inputFile;
    inputFile.open(m_inputFileName, ios::in);

    while (inputFile.fail()){
        cout << "\n\t\t  The input file open failed. Please Check!" << endl;
        cout << "\t>>> Please input the filename of your input again:>>>";
        cin >> m_inputFileName;

        if ((string::npos != m_inputFileName.find("EXIT")) ||
                (string::npos != m_inputFileName.find("QUIT")) ||
                (string::npos != m_inputFileName.find("exit")) ||
                (string::npos != m_inputFileName.find("quit"))) {
            cout << "\t  Program ParallelNonlinearFEM terminated by user." << endl;
            exit(0);
        }

        cout << "\n\t The input filename is " << m_inputFileName << endl;
        inputFile.open(m_inputFileName, ios::in);
    }

    string eachLine = "";
    string inputData = "";
    string eachProp = "";
    string tempLine = "";
    cout << "\n\t Reading data ..." << endl;
    int flag = 0;
    while(getline(inputFile, eachLine)) {
        if ("\r" != eachLine.substr(0, 1)) {
            eachLine.erase(eachLine.find_last_not_of(" ")+1);
            eachLine.erase(0, eachLine.find_first_not_of("\t"));
            eachLine.erase(0, eachLine.find_first_not_of(" "));

            tempLine = eachLine;
            if ("{" == eachLine.substr(0, 1)) {
                flag = flag + 1;
            }
            if ("};" == eachLine.substr(0, 2)) {
                flag = flag - 1;
                if (0 == flag) {
                    eachLine = eachProp.erase(eachProp.find_last_of("\r")) + eachLine;
                }
            }

            if (0 == flag) {
                string temp = eachLine;
                temp.erase(temp.find_last_of("\r"));
                if ('=' != temp[temp.find_last_not_of(" ")]){
                    inputData = inputData + eachLine;
                    lines.emplace_back(eachLine);
                }
                eachLine = tempLine;
                eachProp = eachLine;
            }else{
                eachProp = eachProp.erase(eachProp.find_last_of("\r")) + eachLine;
            }

        }
    }
//    cout << inputData << endl;
    inputFile.close();
}

void InputOutput::ReadTitleInfo(const string &strData)
{
    cout << "\n\t Data File Name:\n\t ==============="  << endl;
    cout << "\t " << strData << endl;
    noTitle = false;
}

void InputOutput::ReadAnalyseTypeInfo(const string &strData)
{
    regex pattern("[1-3]");
    smatch result;
    bool is_smatch = regex_search(strData, result, pattern);
    if (is_smatch) {
        string analTypeStr = result[0];
        sscanf(analTypeStr.c_str(), "%d", &m_analyseType);

        string temp = "";
        for (auto iter = result[0].first; iter != strData.end(); iter++)
        {
            temp = temp + (*iter);
        }
        cout << "\n\t Analysis type ....................................... =  "
             << temp << endl;
        cout << "\t        1 = Plane stress" << endl;
        cout << "\t        2 = Plane strain" << endl;
        cout << "\t        3 = Axisymmetric" << endl;
    }
    else
    {
        cout << "\t ERROR: The analysis Type is Wrong. Please CHECK!" << endl;
        exit(-1);
    }

    if (3 == m_analyseType) {

    }
}

void InputOutput::ReadAxisOfSymmetryInfo(const string &strData)
{
    string temp = "";
    if (string::npos != strData.find("AXIS_OF_SYMMETRY")) {
//        regex pattern("=\s*[XYxy]");
//        smatch result;
//        bool is_smatch = regex_search(strData, result, pattern);
        temp = strData.substr(strData.find("AXIS_OF_SYMMETRY")+16);
        temp.erase(0, temp.find_first_not_of(" "));
        if ("Y" == temp) {
            m_axisOfSymmetry = 1;
        }
        if ("X" == temp) {
            m_axisOfSymmetry = 2;
        }
    }
    else {
        cout << "\n\t Please assign AXIS_OF_SYMMETRY" << endl;
        exit(0);
    }

    cout << "\t Axis of symmetry .................................... =  "
         << temp << endl;
}

void InputOutput::ReadLargeStrainInfo(const string &strData)
{
    string temp = "";
    regex pattern("=[\\s]{0,}[Oo][FNfn][F\\s;][\\s]{0,}");
    smatch result;
    bool is_smatch = regex_search(strData, result, pattern);
    if (is_smatch) {
         temp = result[0];
        if (string::npos != temp.find("O")){
            temp.erase(0, temp.find("O"));
        } else if (string::npos != temp.find("o")) {
            temp.erase(0, temp.find("o"));
        }
        temp.erase(temp.find(";"));
        regex OFF_pattern("[Oo][Ff][Ff]");
        smatch OFF_result;
        bool is_OFF_smatch = regex_match(temp, OFF_result, OFF_pattern);
        if (is_OFF_smatch) {
            m_largeStrain = 0;
        }else {
            m_largeStrain = 1;
        }
    } else {
        temp = "OFF";
        cout << "\t WARNING: Please CHeck the ON-OFF about LARGE_STRAIN_FORMULATION!!" << endl;
    }
    transform(temp.begin(), temp.end(), temp.begin(), ::toupper);
    cout << "\n\t Large deformation flag .............................. =  " << temp << endl;

}

void InputOutput::ReadSolAlgorithmInfo(const string &strData)
{
    regex   pattern("=[\\s]{0,}-{0,1}[1-7]");
    smatch  result{};
    bool    is_match = regex_search(strData, result, pattern);
    if (is_match) {
        string temp = result[0];
        temp.erase(0, 1);
        temp.erase(0, temp.find_first_not_of(" "));
        sscanf(temp.c_str(), "%d", &m_solAlgorithmType);
        cout << "\n\t Nonlinear solution algorithm ........................ =  " << temp << endl;
        cout << "\t        Negative for the arc length method" << endl;
        cout << "\t        1 = Initial stiffness method" << endl;
        cout << "\t        2 = Newton-Raphson tangential stiffness method" << endl;
        cout << "\t        3 = Modified Newton KT1" << endl;
        cout << "\t        4 = Modified Newton KT2" << endl;
        cout << "\t        5 = Secant Newton - Initial stiffness" << endl;
        cout << "\t        6 = Secant Newton - KT1" << endl;
        cout << "\t        7 = Secant Newton - KT2" << endl;
        noSolAlgorithmType = false;
    } else {
        cout << "\t Please Assign Solution Algorithm!!!" << endl;
        exit(0);
    }
}

void InputOutput::ReadArcLengthInfo(const string &strData)
{
    regex   pattern("=[\\s]{0,}[A-Za-Z]{6,9}_[a-zA-Z]\\d{4}");
    smatch  result{};
    bool    is_match = regex_search(strData, result, pattern);
    if (is_match) {
        string temp = result[0];
        temp.erase(0, 1);
        temp.erase(0, temp.find_first_not_of(" "));
        if ("STIFFNESS_SIGN" == temp) {
            m_solArcLengthType = 1;
        }
        if ("SECANT_PATH" == temp) {
            m_solArcLengthType = 2;
        }
    }
    noSolArcLengthType = false;
    cout << "\t Arc-length option ................................... =" << m_solArcLengthType << endl;
    cout << "\t        1 = Follow stiffness determinant sign" << endl;
    cout << "\t        2 = Follow current path" << endl;
}

void InputOutput::ReadElemGroupInfo(const string &strData)
{
    RegexPattern(m_elemGroup, strData);
    m_numElemGroup = m_elemGroup.size();
    cout << "\n\t Element Groups:        Number of element groups =  " << m_numElemGroup << endl;
    cout << "\t ===============" << endl;
    cout << "\t GroupID     Element type    Material type" << endl;
    for (int i = 1; i <= m_numElemGroup; i++) {
        vector<int>&i_elemGroup = m_elemGroup.at(i);
        cout << "\t   " << i_elemGroup[0] << "\t\t  " << i_elemGroup[1] << "\t\t   " << i_elemGroup[2] << endl;
    }
}

void InputOutput::ReadElemTypeInfo(const string &strData)
{
    RegexPattern(m_elemType, strData);

    m_numElemType = m_elemType.size();
    cout << "\n\t Element types:          Number of element types =  " << m_numElemType << endl;
    cout << "\t ==============" << endl;
    cout << "\t ElemTypeID     Elem. Name    Num. Gauss Point" << endl;

    for (int i = 1; i <= m_numElemType; i++) {
        vector<int>&i_elemType = m_elemType.at(i);
        cout << "\t   " << i_elemType[0];
        if (28 == i_elemType[1]){
            cout << "\t\t  QUAD_8";
        }
        cout << "\t\t   " << i_elemType[2] << endl;
    }
}

void InputOutput::ReadGeometryInfo(const string &inputFile)
{
    string filePathName = inputFile;
    filePathName.erase(0, filePathName.find("=")+1);
    filePathName.erase(filePathName.find(";"));
    filePathName.erase(0, filePathName.find_first_not_of(" "));
    FILE *geometryInfo = fopen(filePathName.c_str(), "r");

    regex    space("[\\s]+");
    if (geometryInfo == NULL) {
        cout << "OPEN FILE ERROR!!" << endl;
        exit(0);
    }

    // read node coordinate
    char tempChar[MAX_LENGTH_CHAR] = {};
    while (fgets(tempChar, MAX_LENGTH_CHAR, geometryInfo) != NULL) {
        string temp = tempChar;
        temp.erase(0, temp.find_first_not_of(" "));
        if (string::npos != temp.find("<Nodes>")) {
            temp.erase(0, temp.find_first_not_of(" "));
            if (string::npos != temp.find("CYLINDRICAL"))
            {
                while(fgets(tempChar, MAX_LENGTH_CHAR, geometryInfo) != NULL) {
                    temp = tempChar;
                    if (string::npos != temp.find("</Nodes>")) {
                        break;
                    }
                    int    icount = 0;
                    double floatValue[5] = {};
                    temp.erase(0, temp.find_first_not_of(" "));
                    for (sregex_token_iterator it(temp.begin(), temp.end(), space, -1), end; it != end; it++, icount++) {
                        string coord = it->str();
                        sscanf(coord.c_str(), "%lf", &floatValue[icount]);
                    }
                    int    nodeID  = floatValue[0];
                    double centerX = floatValue[1];
                    double centerY = floatValue[2];
                    double radius  = floatValue[3];
                    double thea    = floatValue[4];
                    double coordX  = centerX + radius * cos(thea);
                    double coordY  = centerY + radius * sin(thea);
                    vector<double> iNodeCoord = {coordX, coordY};
                    m_nodeSet.insert(pair<int, vector<double>>(nodeID, iNodeCoord));
                }
            } else if (string::npos != temp.find("CARTESIAN"))
            {
                while (fgets(tempChar, MAX_LENGTH_CHAR, geometryInfo)) {
                    temp = tempChar;
                    if (string::npos != temp.find("</Nodes>")) {
                        break;
                    }
                    int icount = 0;
                    double floatValue[3];
                    for (sregex_token_iterator it(temp.begin(), temp.end(), space, -1), end;
                         it != end; it++, icount++) {
                        cout << it->str() << endl;
                        string coord = it->str();
                        sscanf(coord.c_str(), "%lf", &floatValue[icount]);
                    }
                    int nodeID = floatValue[0];
                    vector<double> iNodeCoord = {floatValue[1], floatValue[2]};
                    m_nodeSet.insert(pair<int, vector<double>>(nodeID, iNodeCoord));
                }
            }
            m_numNode = m_nodeSet.size();
            cout << "\n\t Nodal point co-ordinates:       Number of nodes =  " << m_numNode << endl;

        }
        // read element connection
        else if (string::npos != temp.find("<Elements>")) {
            while (fgets(tempChar, MAX_LENGTH_CHAR, geometryInfo)) {
                temp = tempChar;
                if (string::npos != temp.find("</Elements>")) {
                    break;
                }

                vector<int> ielemConnect = {};
                temp.erase(0, temp.find_first_not_of(" "));
                for (sregex_token_iterator it(temp.begin(), temp.end(), space, -1), end; it != end; it++) {
                    string coord = it->str();
                    int intValue = 0;
                    sscanf(coord.c_str(), "%d", &intValue);
                    ielemConnect.emplace_back(intValue);
                }
                m_elemConnect.insert(pair<int, vector<int>>(ielemConnect[0], ielemConnect));
                m_elemConnect.at(ielemConnect[0]).erase(m_elemConnect.at(ielemConnect[0]).begin());
            }
            m_numElem = m_elemConnect.size();
            cout << "\n\t Element connectivities:      Number of elements =  " << m_numElem << endl;
        }
        //read boundary condition
        else if (string::npos != temp.find("<NodeDispBoundary>")) {
            while (fgets(tempChar, MAX_LENGTH_CHAR, geometryInfo)) {
                temp = tempChar;
                if (string::npos != temp.find("</NodeDispBoundary>")) {
                    break;
                }
                vector<string> inodeDispStr = {};
                temp.erase(0, temp.find_first_not_of(" "));
                for (sregex_token_iterator it(temp.begin(), temp.end(), space, -1), end; it != end; it++) {
                    string tempStr = it->str();
                    inodeDispStr.emplace_back(tempStr);
                }
                int nodeID = 0;
                sscanf(inodeDispStr[0].c_str(), "%d", &nodeID);
                nodeBoundaryCondition inodeDisp = {};
                inodeDisp.Dof = inodeDispStr[1];
                for (int i = 2; i < inodeDispStr.size(); i++) {
                    double doubleValue = 0.;
                    sscanf(inodeDispStr[i].c_str(), "%lf", &doubleValue);
                    inodeDisp.value.emplace_back(doubleValue);
                }
                m_prescribeNodeDisp.insert(pair<int, nodeBoundaryCondition>(nodeID, inodeDisp));
            }
        }
        // read node force
        else if (string::npos != temp.find("<NodeForceBoundary>")) {
            while (fgets(tempChar, MAX_LENGTH_CHAR, geometryInfo)) {
                temp = tempChar;
                if (string::npos != temp.find("</NodeForceBoundary>")) {
                    break;
                }
                vector<string> inodeForceStr;
                for (sregex_token_iterator it(temp.begin(), temp.end(), space, -1), end; it != end; it++) {
                    string tempStr = it->str();
                    inodeForceStr.emplace_back(tempStr);
                }
                int nodeID = 0;
                sscanf(inodeForceStr[0].c_str(), "%d", &nodeID);
                nodeBoundaryCondition inodeForce = {};
                inodeForce.Dof = inodeForceStr[1];
                for (int i = 2; i < inodeForceStr.size(); i++) {
                    double doubleValue = 0.;
                    sscanf(inodeForceStr[i].c_str(), "%lf", &doubleValue);
                    inodeForce.value.emplace_back(doubleValue);
                }
                m_prescribeNodeForce.insert(pair<int, nodeBoundaryCondition>(nodeID, inodeForce));
            }
        }
        // read distribute force
        else if (string::npos != temp.find("<DistributeForce>")) {
            while (fgets(tempChar, MAX_LENGTH_CHAR, geometryInfo)) {
                temp =tempChar;
                if (string::npos != temp.find("</DistributeForce>")) {
                    break;
                }
                vector<string> inodeDistributeForceStr;
                for (sregex_token_iterator it(temp.begin(), temp.end(), space, -1), end; it != end; it++) {
                    string tempStr = it->str();
                    inodeDistributeForceStr.emplace_back(tempStr);
                }

                int elemID = 0;
                sscanf(inodeDistributeForceStr[0].c_str(), "%d", &elemID);
                elemDistributeForce inodeDistributeForce = {};
                sscanf(inodeDistributeForceStr[1].c_str(), "%d", &inodeDistributeForce.nodeID);
                for (int i = 2; i < inodeDistributeForceStr.size(); i++) {
                    double doubleValue = 0.;
                    sscanf(inodeDistributeForceStr[i].c_str(), "%lf", &doubleValue);
                    inodeDistributeForce.value.emplace_back(doubleValue);
                }
                if (0 == m_distributeForce.count(elemID))
                {
                    m_distributeForce.insert(pair<int, vector<elemDistributeForce>>(elemID, {}));
                }
                vector<elemDistributeForce> &ielemDistributeForce = m_distributeForce.at(elemID);
                ielemDistributeForce.emplace_back(inodeDistributeForce);
            }
        }

    }
    fclose(geometryInfo);
}

void InputOutput::ReadNewtonRaphsonAlgorithmInfo(const string &strData)
{
    //cout << strData << endl;
    regex   pattern("=[\\s]*\\[(.*)\\];");
    smatch  result ={};


    if (regex_search(strData, result, pattern)) {
        string temp = result[0];
        temp.erase(0, temp.find_first_of(" "));
        temp.erase(0, temp.find_first_not_of(" "));
        temp.erase(temp.size()-2);
        temp.erase(0, 1);

        regex   semicolon(";");
        regex   space("[\\s]+");
        for (sregex_token_iterator it(temp.begin(), temp.end(), semicolon, -1), end; it != end; it++) {
            string iIncrement = it->str();
            vector<double> i_increments = {};
            int iicount = 0;
            for (sregex_token_iterator iit(iIncrement.begin(), iIncrement.end(), space, -1), end; iit != end; iit++, iicount++) {
                string tempStr = iit->str();
                double doubleValue = 0.;
                sscanf(tempStr.c_str(), "%lf", &doubleValue);
                i_increments.emplace_back(doubleValue);
            }
            m_increments.emplace_back(i_increments);
        }

    }
}

void InputOutput::ReadDispBoundaryInfo(const string &strData)
{

}

void InputOutput::ReadMaterialInfo(const string &strData)
{
    regex   pattern("[\\w]*[1-9][0-9]*[\\s]*=[\\s]*\\{(.*)};");
    smatch  result = {};
    string::const_iterator iter_begin = strData.cbegin();
    string::const_iterator iter_end   = strData.cend();
    while (regex_search(iter_begin, iter_end, result, pattern)) {
//        cout << result[0] << endl;
        string material = result[0];
        regex  temp_pattern("[tT][yY][pP][eE][\\s]*=[\\s]*[\\w]*;");
        smatch temp_result   = {};
        bool   is_temp_match = regex_search(material, temp_result, temp_pattern);
        if (is_temp_match) {
            string mat_type = temp_result[0];
            mat_type.erase(0, mat_type.find("=")+1);
            mat_type.erase(mat_type.find(";"));
            mat_type.erase(0, mat_type.find_first_not_of(" "));

            ReadMaterialPara(mat_type, material);
        } else {
            cout << "Please Assign Material Type!!!" << endl;
            exit(0);
        }
        iter_begin = result[0].second;
    }
}

void InputOutput::ReadMaterialPara(const string &matType, const string &para)
{
    if ("VON_MISES" == matType) {

    } else if (true) {

    }
}

void InputOutput::CheckInputFile()
{
    if (noTitle) {
        cout << "\t WARNING: The inpu file not include TITLE. Please check!" << endl;

    } else if ((1 != m_analyseType) && (2 != m_analyseType) && (3 != m_analyseType)) {
        cout << "\t ERROR: The analysis Type is Wrong. Please CHECK!" << endl;
        exit(0);
    } else if (noSolAlgorithmType) {
        cout << "\t ERROR: Please assign the SOL ALGORITHM in input file!" << endl;
        exit(0);
    } else if (m_numElemGroup <= 0) {
        cout << "\t" << endl;
        exit(0);
    }

}

void InputOutput::RegexPattern(unordered_map<int, vector<int> > &tempMap, const string &strData)
{
    regex   pattern("[\\w]*[1-9][0-9]*[\\s]*=\\{([\\s]*[\\w]*[\\s]*=[\\s]*[1-9][0-9]*;)*\\};");
    smatch  result = {};
    string::const_iterator iter_begin = strData.cbegin();
    string::const_iterator iter_end   = strData.cend();
    while (regex_search(iter_begin, iter_end, result, pattern)) {
        string  temp = result[0];
        regex   temp_pattern("=[\\s]*[1-9][0-9]*");
        smatch  temp_result{};
        string::const_iterator temp_iter_begin = temp.cbegin();
        string::const_iterator temp_iter_end   = temp.cend();
        vector<int> i_vector{};
        while(regex_search(temp_iter_begin, temp_iter_end, temp_result, temp_pattern)) {
            string temp_str = temp_result[0];
            temp_str.erase(0, 1);
            temp_str.erase(0, temp_str.find_first_not_of(" "));
            int i_value = 0;
            sscanf(temp_str.c_str(), "%d", &i_value);
            i_vector.emplace_back(i_value);
            temp_iter_begin = temp_result[0].second;
        }
        tempMap.insert(pair<int, vector<int>>(i_vector[0], i_vector));
        iter_begin = result[0].second;
    }
}


int InputOutput::GetNodeNum()
{
    m_numNode = m_nodeSet.size();
    return m_numNode;
}

int InputOutput::GetElemNum()
{
    m_numElem = m_elemConnect.size();
    return m_numElem;
}

int InputOutput::GetAnalyseType()
{
    return m_analyseType;
}

int InputOutput::GetAxisOfSymmetry()
{
    return m_axisOfSymmetry;
}

int InputOutput::GetLargeStrain()
{
    return m_largeStrain;
}

int InputOutput::GetSolAlgorithmType()
{
    return m_solAlgorithmType;
}

int InputOutput::GetSolArcLengthType()
{
    return m_solArcLengthType;
}

unordered_map<int, vector<double>> InputOutput::GetNodeCoord()
{
    return m_nodeSet;
}

unordered_map<int, vector<int>> InputOutput::GetElemType()
{
    return m_elemType;
}

unordered_map<int, vector<int>> InputOutput::GetElemGroup()
{
    return m_elemGroup;
}

unordered_map<int, nodeBoundaryCondition> InputOutput::GetPrescribeNodeDisp()
{
    return m_prescribeNodeDisp;
}

unordered_map<int, nodeBoundaryCondition> InputOutput::GetPrescribeNodeForce()
{
    return m_prescribeNodeForce;
}

unordered_map<int, vector<int>> InputOutput::GetElemConnect()
{
    return m_elemConnect;
}

unordered_map<int, vector<elemDistributeForce>> InputOutput::GetDistributeForce()
{
    return m_distributeForce;
}

vector<vector<double>> InputOutput::GetIncrements()
{
    return m_increments;
}
