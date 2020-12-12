#include "Mesher.h"
#include "LammpsInterface.h"
#include "Utility.h"
#include "Function.h"
#include "MPI_Wrappers.h"

using namespace std;
using ATC_Utility::is_numeric;
using ATC_Utility::to_string;;
using ATC_Utility::parse_direction;;
using MPI_Wrappers::print_msg;
using MPI_Wrappers::print_msg_once;

const double _pi = 4.0*atan(1.0);

namespace ATC {
  Mesher::Mesher():
    mesh_(nullptr),
    nNodes_(0),
    nElements_(0),
    nNodeSets_(0),
    nSideSets_(0),
    nElemSets_(0)
  {
  };


  //-----------------------------------------------------------------
  bool Mesher::modify(int narg, char **arg)
  //-----------------------------------------------------------------
  {
    bool match = false;
    /*! \page man_mesh_create fix_modify AtC mesh create
      \section syntax
      fix_modify AtC mesh create <nx> <ny> <nz> <region-id> 
      <f|p> <f|p> <f|p> \n
      - nx ny nz = number of elements in x, y, z
      - region-id = id of region that is to be meshed
      - f p p = periodicity flags for x, y, z
      fix_modify AtC mesh create <nx> <ny> <nz> <region-id> dx <function_name parameters|s_1 s_2 ... s_nx|position-number-density list>
      fix_modify AtC mesh create prism <nvertices> <R> <axis x|y|z> <zlo> <zhi> 
      <nside> <nheight> [ <x0> <y0> ] \n
      \section examples
      <TT> fix_modify AtC mesh create 10 1  1  feRegion p p p </TT> \n
      <TT> fix_modify ATC mesh create ${nx} ${ny} 1  BOX f f p dy position-number-density 0 0 4.0 ${w1} ${n1} 1. ${w2} ${n2} 1. 1 ${ny} 4.0 </TT> \n
      \section description
      Creates a uniform|rectangular|prismatic mesh in a region
      \section restrictions
      Creates only uniform rectangular or graded rectangular grids in a rectangular region or prismatic mesh of hexahedral elements in prismatic region
      \section related
      \ref man_mesh_quadrature
      \section default
      When created, mesh defaults to gauss2 (2-point Gaussian) quadrature. 
      Use "mesh quadrature" command to change quadrature style.
    */
    int argIdx = 0;
    if (strcmp(arg[argIdx],"prism")==0) {
      int nvertices = atoi(arg[++argIdx]);
      double radius = atof(arg[++argIdx]);
      int axis = parse_direction(arg[++argIdx]);
      double zlo = atof(arg[++argIdx]);
      double zhi = atof(arg[++argIdx]);
      int nside = atoi(arg[++argIdx]);
      int nheight = atoi(arg[++argIdx]);
      bool zperiodicity = false;
      if (argIdx+1 < narg) 
        { zperiodicity = (strcmp(arg[++argIdx],"p")==0) ? true : false;}
      double x0 = 0.;
      double y0 = 0.;
      if ( narg-1 > argIdx) {
        x0 = atof(arg[++argIdx]);
        y0 = atof(arg[++argIdx]);
      }
      create_prismatic_mesh(nvertices, radius, axis, zlo, zhi, nside, nside, nheight, zperiodicity, x0, y0); 
      match = true;
    }
    else { // rectangular
      int nx = atoi(arg[argIdx++]);
      int ny = atoi(arg[argIdx++]);
      int nz = atoi(arg[argIdx++]);
      string box = arg[argIdx++];

      Array<bool> periodicity(3);
      periodicity(0) = (strcmp(arg[argIdx++],"p")==0) ? true : false; 
      periodicity(1) = (strcmp(arg[argIdx++],"p")==0) ? true : false;
      periodicity(2) = (strcmp(arg[argIdx++],"p")==0) ? true : false;

      if (argIdx < narg ) {
        create_rectangular_mesh(nx, ny, nz, box.c_str(), periodicity, 
                                argIdx,narg,arg);
        match = true;
      }
      else  {
        create_uniform_mesh(nx, ny, nz, box.c_str(), periodicity);
        match = true;
      }
    }
    return match;
  }
  //-----------------------------------------------------------------
  // create a uniform rectangular mesh on a rectangular region
  //-----------------------------------------------------------------
  void Mesher::create_uniform_mesh(int nx, int ny, int nz, 
     const char * regionName, Array<bool> periodicity)
  {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xscale, yscale, zscale;

    // check to see if region exists and is it a box, and if so get the bounds
    bool isBox;
    isBox = ATC::LammpsInterface::instance()->region_bounds(
                                    regionName,
                                    xmin, xmax,
                                    ymin, ymax,
                                    zmin, zmax,
                                    xscale, yscale, zscale);
    if (!isBox) throw ATC_Error("Region for FE mesh is not a box");
    
    mesh_ = new FE_Uniform3DMesh(nx,ny,nz,
                                   xmin, xmax,
                                   ymin, ymax,
                                   zmin, zmax,
                                   periodicity,
                                   xscale, yscale, zscale);
  }

  //-----------------------------------------------------------------
  // create a non-uniform rectangular mesh on a rectangular region
  //-----------------------------------------------------------------
  void Mesher::create_rectangular_mesh(int nx, int ny, int nz,
    const char * regionName,
    Array<bool> periodicity,
    int argIdx, int narg, char **arg)
  {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xscale, yscale, zscale;

    // check to see if region exists and is it a box, and if so get the bounds
    bool isBox;
    isBox = ATC::LammpsInterface::instance()->region_bounds(
                                    regionName,
                                    xmin, xmax,
                                    ymin, ymax,
                                    zmin, zmax,
                                    xscale,
                                    yscale,
                                    zscale);
    if (!isBox) throw ATC_Error("Region for FE mesh is not a box");

    Array<double> dx(nx),dy(ny),dz(nz);
    dx = 0;
    dy = 0;
    dz = 0;
    if (dx(0) == 0.0) dx = (xmax-xmin)/dx.size();
    if (dy(0) == 0.0) dy = (ymax-ymin)/dy.size();
    if (dz(0) == 0.0) dz = (zmax-zmin)/dz.size();
    while (argIdx < narg) {
      if      (strcmp(arg[argIdx],"dx")==0) {
        parse_partitions(++argIdx, narg, arg, 0, dx, xmin, xmax);
      }
      else if (strcmp(arg[argIdx],"dy")==0) {
        parse_partitions(++argIdx, narg, arg, 1, dy, ymin, ymax);
      }
      else if (strcmp(arg[argIdx],"dz")==0) {
        parse_partitions(++argIdx, narg, arg, 2, dz, zmin, zmax);
      }
    }

    mesh_ = new FE_Rectangular3DMesh(dx,dy,dz,
                                       xmin, xmax,
                                       ymin, ymax,
                                       zmin, zmax,
                                       periodicity,
                                       xscale, yscale, zscale);
    stringstream ss;
#ifdef ATC_VERBOSE
    print_partitions(xmin,xmax,dx);
    print_partitions(ymin,ymax,dy);
    print_partitions(zmin,zmax,dz);
#endif
  }
  //-----------------------------------------------------------------
  void Mesher::parse_partitions(int & argIdx, int narg, char ** arg, 
    int idof, Array<double> & dx, double xmin, double xmax ) const
  //-----------------------------------------------------------------
  {
    double x[3] = {0,0,0}; 
    int nx = dx.size();
    // parse relative values for each element
    if (is_numeric(arg[argIdx])) {
      for (int i = 0; i < nx; ++i) { 
        if (is_numeric(arg[argIdx])) { dx(i) = atof(arg[argIdx++]); }
        else throw ATC_Error("not enough element partitions");
      }
    }
    // each segment of the piecewise funcion is length-normalized separately
    else if (strcmp(arg[argIdx],"position-number-density")==0) { 
      argIdx++;
      double y[nx],w[nx];
      int n[nx];
      int nn = 0;
      while (argIdx < narg) { 
        if (! is_numeric(arg[argIdx])) break;
        y[nn]   = atof(arg[argIdx++]);
        n[nn]   = atoi(arg[argIdx++]);
        w[nn++] = atof(arg[argIdx++]);
      }
      if (n[nn-1] != nx)  throw ATC_Error("total element partitions do not match");
      int k = 0;
      for (int i = 1; i < nn; ++i) { 
        int dn = n[i]-n[i-1];
        double dy = y[i]-y[i-1];
        double w0 = w[i-1];
        double dw = w[i]-w0;
        double lx = 0;
        double l[dn];
        for (int j = 0; j < dn; ++j) {
          double x = (j+0.5)/dn; 
          double dl = w0+x*dw;
          lx += dl;
          l[j]= dl;
        }
        double scale = dy/lx;
        for (int j = 0; j < dn; ++j) {
          dx(k++) = scale*l[j];
        }
      }
    }
    // construct relative values from a density function
    // evaluate for a domain (0,1)
    else {
      XT_Function * f = XT_Function_Mgr::instance()->function(&(arg[argIdx]),narg-argIdx);
      argIdx++;
      while (argIdx < narg) { 
        if (! is_numeric(arg[argIdx])) break;
        argIdx++;
      }
      double diffx = xmax-xmin;
      for (int i = 0; i < nx; ++i) { 
        x[idof] = xmin+diffx*(i+0.5)/nx; dx(i) = f->f(x,0.); 
      }
    }
    for (int i = 0; i < nx; ++i) { cout << dx(i) << " "; }; cout << "\n";
  }
  //-----------------------------------------------------------------
  // 
  //-----------------------------------------------------------------
  void Mesher::print_partitions(double xmin, double xmax, 
    Array<double> & dx) const
  {
    stringstream msg;
    msg.precision(3);
    msg << std::fixed;
    msg <<  "\nindex weight fraction location size[A] size[uc]:\n";
    double sum = 0;
    for (int i = 0; i < dx.size(); ++i) { sum += dx(i); }
    double xn= 1.0/sum;
    double xs= xn*(xmax-xmin);
    double xl= xs/(ATC::LammpsInterface::instance()->xlattice());
    double sumn = 0, sums = 0, suml = 0;
    double x = xmin;
    for (int i = 0; i < dx.size(); ++i) { 
       double dxn = dx(i)*xn; sumn += dxn;
       double dxs = dx(i)*xs; sums += dxs;
       double dxl = dx(i)*xl; suml += dxl;
       msg << std::setw(5) << i+1
           << std::setw(7) << dx(i)
           << std::setw(9) << dxn
           << std::setw(9) << x
           << std::setw(8) << dxs
           << std::setw(9) << dxl << "\n";
       x += dxs;
     }
     msg << "sum  " << setw(7) << sum
                    << setw(9) << sumn
                    << setw(9) << x
                    << setw(8) << sums
                    << setw(9) << suml << "\n";

  }

  //------------------------------------------------------------------
  // create a radial mesh in a right polygonal prism
  //------------------------------------------------------------------
  void Mesher::create_prismatic_mesh(
    int nvertices,
    double radius,
    int axis,
    double zlo, double zhi,
    int nradial, int nside, int nheight,
    bool zperiodicity,
    double center1, double center2)
  {
    if (axis != 2) throw ATC_Error("only z-aligned prismatic mesh supported");
    double R = radius;
    vector< DENS_VEC > vertices;
    for (int i=0; i < nvertices; ++i) {
      double theta = 2.*_pi*(i+0.5)/nvertices;
      DENS_VEC v(3);
      v(0) = R*sin(theta);
      v(1) = R*cos(theta);
      v(2) = 0;
      vertices.push_back(v);
    }
    nNodes_ = nvertices*(nside+1)*(nside+1)*(nheight+1);
    nElements_ = nvertices*nside*nside*nheight;

    // nodes
    nodeCoords_.reset(3, nNodes_, false);
    double h = zlo;
    double s = 1./nside;
    DENS_VEC c(3);
    c[0] = center1; c[1] = center2; c[2] = h;
    DENS_VEC p0(3),p1(3),p2(3),p3(3);
    DENS_VEC x0(3),x1(3),d1(3),d2(3);
    DENS_VEC x(3), dx(3);
    int n = 0;
    for (int I=0; I < nvertices; ++I) {
      int J = I - 1; if (J < 0) J = nvertices-1;
      int K = I + 1; if (K > nvertices-1) K = 0;
      DENS_VEC & v0 = vertices[J];
      DENS_VEC & v1 = vertices[I];
      DENS_VEC & v2 = vertices[K];
      DENS_VEC & v3 = c;
      p1 = v1;               p1(2) = h;
      p0 = v1 + 0.5*(v0-v1); p0(2) = h;
      p2 = v1 + 0.5*(v2-v1); p2(2) = h;
      p3 = v3;               p3(2) = h;
      d1 = s*(p0-p1);
      d2 = s*(p3-p2);
      for (int i = 0; i < nside; ++i) { 
        double ii = i;
        x0 = p1+ii*d1;
        x1 = p2+ii*d2;
        dx = s*(x1-x0);
        for (int j = 0; j < nside+1; ++j){
          double jj = j;
          x = x0+jj*dx; 
          nodeCoords_(0,n) = x(0);
          nodeCoords_(1,n) = x(1);
          nodeCoords_(2,n) = h;
        }
        n++;
      }
    }
    nodeCoords_(0,n) = c(0);
    nodeCoords_(1,n) = c(1);
    nodeCoords_(2,n) = c(2);
    n++;
    int p = n;
    double dz = (zhi-zlo)/nheight;
    for (int k = 0; k < nheight; ++k) { 
      h += dz;
      for (int i = 0; i < p; ++i) {
        nodeCoords_(0,n) = nodeCoords_(0,i);
        nodeCoords_(1,n) = nodeCoords_(1,i);
        nodeCoords_(2,n) = h;
        n++;
      }
    }

    // elements
    conn_.reset(8, nElements_); // HEX8 only
    n = 0;
    int ss = nside*(nside+1);
    for (int k = 0; k < nside; ++k) { 
      for (int I = 0;  I < nvertices; ++I) {
        int ii = I*ss+k*p;
        for (int j = 0; j < nside-1; ++j) { 
          for (int i = 0;  i < nside ; ++i ) {
            ii++;
            int jj = ii + nside + 1;
            int kk = ii + p;
            int ll = jj + p;
            conn_(n,0) = ii;
            conn_(n,1) = jj;
            conn_(n,2) = jj+1;
            conn_(n,3) = ii+1;
            conn_(n,4) = kk;
            conn_(n,5) = ll;
            conn_(n,6) = ll+1;
            conn_(n,7) = kk+1;
            n++;
          }
          ii++;
        }
        int jj = (I-1)*ss+k*p;
        if (I == 0) jj = (n-1)*ss+k*p;
        for (int i = 0; i < nside; ++i) { 
          ii++;
          int i2 = ii + 1;
          jj += nside+1;
          int j2 = jj+nside+1;
          if (i == nside-1) j2 = (k+1)*p;
          conn_(n,0) = ii;
          conn_(n,1) = jj;
          conn_(n,2) = j2;
          conn_(n,3) = i2;
          conn_(n,4) = ii+p;
          conn_(n,5) = jj+p;
          conn_(n,6) = j2+p;
          conn_(n,7) = i2+p;
        }
      }
    }
    nodeSets_.reset(1);
    nodeSets_(0).first = "all";
    for (int i = 0; i < nNodes_; ++i) {
      nodeSets_(0).second.insert(i);
    }

    Array<bool> periodicity(3);
    periodicity(0) = false;
    periodicity(1) = false;
    periodicity(2) = zperiodicity;
    mesh_  = new FE_3DMesh("HEX8",
                         nNodes_, nElements_,
                         &conn_, &nodeCoords_,
                         periodicity,
                         &nodeSets_);
  };
}
