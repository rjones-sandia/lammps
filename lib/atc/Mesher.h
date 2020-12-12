#ifndef MESHER_H
#define MESHER_H

#include "Array2D.h"
#include "MatrixDef.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"
#include "FE_Mesh.h"

namespace ATC {

  class Mesher {

  public:
    Mesher();
    ~Mesher(){};
    bool modify(int narg, char **arg);
    FE_Mesh * mesh() { return mesh_; }
    
  private:
    /** create a uniform structured mesh */
    void create_uniform_mesh(
                     int nx, 
                     int ny, 
                     int nz,
                     const char *regionName,
                     Array<bool> periodic);
    /** create a rectangular structured mesh */
    void create_rectangular_mesh(
                     int nx, 
                     int ny, 
                     int nz,
                     const char *regionName,
                     Array<bool> periodic,
                     int iarg, int narg, char **args);


    /** create a primatic mesh */
    void create_prismatic_mesh(
                     int nvertices,
                     double radius,
                     int axis,
                     double zlo, double zhi,
                     int nradial, int nside, int nheight,
                     bool zperiodicity,
                     double center1, double center2);

    /** helpers */
    void parse_partitions(int & argIdx, int narg, char ** arg,
      int idof, Array<double> & dx, double xmin, double xmax ) const;
    void print_partitions(double xmin, double xmax, Array<double> &dx) const;

 
    /** data */
    FE_Mesh * mesh_;
    int nNodes_;
    int nElements_;
    int nNodeSets_;
    int nSideSets_;
    int nElemSets_;
    ATC_matrix::Array2D<int> conn_;
    DENS_MAT  nodeCoords_;
    ATC_matrix::Array<std::pair<std::string,std::set<int> > >  nodeSets_;
    ATC_matrix::Array<std::pair<std::string,std::set<std::pair<int,int> > > > sideSets_;
    ATC_matrix::Array<std::pair<std::string,std::set<int> > > elemSets_;
  };

};
#endif
