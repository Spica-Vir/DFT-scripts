#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <map>
#include <vector>
#include <math.h>

using namespace std;

class CubeFile
{
    public:
        /* 
        npt: number of data points along abc
        latt: lattice in AU
        step: Vector of a single step along abc
        ori: Origin, in AU
        data: A value at a grid point only
        */
        const std::string cubeFileName;
        std::string title;
        int nAtoms, npt[3]; 
        double latt[3][3], step[3][3], ori[3];
        std::vector<std::vector<double>> atomCoords;
        std::vector<int> atomSymbols;
        std::vector<float> data;

        bool readCUBE_3D(const char* cubeFileName);
        bool writeXSF_3D(const char* xsfFileName, const char* datatype);

    private:
        // constants
        float const br2ang=0.5291772083; // Par(32)
        float const au2v=27.2113834; // Par(3)

        // data conversion types (seems to have error with map)
        const string chgtype="charge", pottype="potential", noconvert="unknown";
};

// Convert from CUBE to XSF
int CUBEtoXSF_3D (const char* cubeFileName, const char* xsfFileName, const char* datatype);
