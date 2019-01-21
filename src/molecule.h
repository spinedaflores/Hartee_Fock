#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
 
class Molecule
{
  public:
    int natom;
    int charge;
    int* zvals;
    double** geom;
    std::string point_group;
 
    void print_geom_data();

    void rotate(double phi);
    void translate(double x, double y, double z);

    double bond(int atom1, int atom2);
    double angle(int atom1, int atom2, int atom3);
    double oop(int atom_i,int atom_j,int atom_k,int atom_l);
    double torsion(int atom_i,int atom_j,int atom_k,int atom_l);

    double unit(int cart, int i, int j);

 
    Molecule(const char* const filename, int const q);
    ~Molecule();
};

#endif
