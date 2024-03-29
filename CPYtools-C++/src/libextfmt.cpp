#include "libextfmt.h"

/*
H.Z. would like to express his sincere thankfulness to ChatGPT
though most of code is terribly wrong.
*/

bool CubeFile::readCUBE_3D(const char* cubeFileName)
{
    std::ifstream cubeFile(cubeFileName);
    if (!cubeFile.is_open()) {
        return false;
    }

    // title
    std::string line;
    std::getline(cubeFile, title);
    size_t start = title.find_first_not_of(" \t");
    size_t end = title.find_last_not_of(" \t");
    title = title.substr(start, end - start + 1);
    std::getline(cubeFile, line);

    // a periodic box file is read anyway
    cubeFile >> nAtoms >> ori[0] >> ori[1] >> ori[2];
    for (int i = 0; i < 3; ++i) {
        ori[i] = ori[i];
    }

    for (int i = 0; i < 3; ++i) {
        cubeFile >> npt[i] >> step[i][0] >> step[i][1] >> step[i][2];
        for (int j = 0; j < 3; ++j) {
            latt[i][j] = step[i][j] * (npt[i] - 1);
        }
    }

    // Read atomic coordinates
    double charge, cord1, cord2, cord3;
    int symbol;
    for (int i = 0; i < nAtoms; ++i) {
        cubeFile >> symbol >> charge >> cord1 >> cord2 >> cord3;
        atomSymbols.push_back(symbol);
        atomCoords.push_back({cord1, cord2, cord3});
    }

    // Read 3D data grid
    int nx=npt[0], ny=npt[1], nz=npt[2];
    float value;
    while (cubeFile >> value) {
        data.push_back(value);
    }
    cubeFile.close();

    return true;
}

bool CubeFile::writeXSF_3D(const char* xsfFileName, const char* datatype)
/*
datatype: Accepted values
    charge: Bohr^-3 --> Angstrom^-3
    potential: Hartree/e --> V
    none: No convert
*/
{
    std::ofstream xsfFile(xsfFileName);
    if (!xsfFile.is_open()) {
        return false;
    }

    // Write XSF header
    xsfFile << "# Generated by CRYSTALpytools\n\n"
            << "CRYSTAL\n"
            << "PRIMVEC\n";
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            xsfFile << fixed << setprecision(6) << setw(12) << right << latt[i][j]*br2ang;
        }
        xsfFile << "\n";
    }

    // Write atomic coordinates to XSF file
    xsfFile << "PRIMCOORD\n"
            << nAtoms << " 1\n";
    for (int i = 0; i < nAtoms; ++i) {
        xsfFile << setw(5) << left << atomSymbols[i];
        for (int j = 0; j < 3; ++j) {
            xsfFile << fixed << setprecision(6) << setw(12) << right << (atomCoords[i][j] - ori[j])*br2ang; // Shift structure to data grid origin - for vesta
        }
        xsfFile << "\n";
    }

    // Compute spanning vector
    float scale[3];

    for (int i=0; i<3; ++i) {
        scale[i] = pow(pow(step[i][0], 2) + pow(step[i][1], 2) + pow(step[i][2], 2), 0.5) / pow(pow(latt[i][0], 2) + pow(latt[i][1], 2) + pow(latt[i][2], 2), 0.5);
        scale[i] = scale[i] * (npt[i] - 1);
    }

    // Read and write the cube data
    int nx=npt[0], ny=npt[1], nz=npt[2];
    std::string comment;

    comment = title;
    std::replace(title.begin(), title.end(), ' ', '_');
    xsfFile << "\nBEGIN_BLOCK_DATAGRID_3D\n"
            << "  " << title << "\n"
            << "  BEGIN_DATAGRID_3D_UNKNOWN\n"
            << "    " << setw(8) << nx << setw(8) << ny << setw(8) << nz << "\n";
    xsfFile << fixed << setprecision(6) << noshowpos;
    xsfFile << "    "
            << setw(15) << right << 0.
            << setw(15) << right << 0.
            << setw(15) << right << 0. << "\n";
    xsfFile << fixed << setprecision(6) << noshowpos;
    for (int i=0; i<3; ++i) {
        xsfFile << "    ";
        for (int j=0; j<3; ++j) {
            xsfFile << setw(15) << right << latt[i][j] * scale[i] * br2ang;
        }
        xsfFile << "\n";
    }

    title = comment;

    // C order to Fortran order
    int totaldata=nx*ny*nz, idx, counter=0;
    float value, convert;

    if (datatype == chgtype) {
        convert = 1 / pow(br2ang, 3);
    }
    else if (datatype == pottype) {
        convert = au2v;
    }
    else if (datatype == noconvert) {
        convert = 1.0;
    }
    else {
        return false;
    }

    xsfFile << scientific << setprecision(6) << noshowpos;
    for (int k=0; k<nz; ++k){
        for (int j=0; j<ny; ++j){
            for  (int i=0; i<nx; ++i) {
                idx = i * ny * nz + j * nz + k;
                value = data[idx] * convert;
                xsfFile << setw(16) << right << value;
                if (counter%5 == 4) {
                    xsfFile << "\n";
                }
                counter++;
            }
        }
    }
    if (--counter%5 != 4) {
        xsfFile << "\n";
    }
    xsfFile << "  END_DATAGRID_3D_UNKNOWN\n";
    xsfFile << "END_BLOCK_DATAGRID_3D\n";
    xsfFile.close();

    return true;
}

// Convert from CUBE to XSF
int CUBEtoXSF_3D (const char* cubeFileName, const char* xsfFileName, const char* datatype) {
    CubeFile cube;
    bool read, write;

    read = cube.readCUBE_3D(cubeFileName);
    if (read != true) {
        int result = 1;
        return result;
    }

    write = cube.writeXSF_3D(xsfFileName, datatype);
    if (write != true) {
        int result = 2;
        return result;
    }

    int result = 0;
    return result;
}
