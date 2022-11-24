#ifndef MATERIALLIB_H
#define MATERIALLIB_H

class MaterialLib
{
public:
    MaterialLib();
    ~MaterialLib();

public:
    double *m_D;
    double *m_stress;
};

class LinearElasticity : public MaterialLib {

};

class PlainStrain
{
public:
    PlainStrain();
};

#endif // MATERIALLIB_H
