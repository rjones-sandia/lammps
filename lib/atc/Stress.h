#ifndef STRESS_H
#define STRESS_H

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"
#include "ATC_TypeDefs.h"
#include "NonLinearSolver.h"

namespace ATC {
  enum ElasticityTensorType {FIRST_ELASTICITY_TENSOR=0, SECOND_ELASTICITY_TENSOR};
  class StressManager
  {
    public:
      StressManager() {};
      virtual ~StressManager() {};
      class Stress * create(std::string matParamFile);
  };
  /**
   * @class Stress
   * @brief Base class that defines interface for a constitutive law 
   * @brief that computes stress given all field and gradient information.
   */
  class Stress
  {
    public:
      Stress()  {};
      virtual ~Stress() {};
      virtual void initialize(void){};
      //* Returns parameter values, (Nothing uses this).
      virtual void parameters(std::map<std::string,double> &parameters) {}
      //* Computes stress given a displacement gradient.
      //* Units: mvv/L^3 (i.e. for units Real: g/(mol ps^2 A^2) )
      virtual void stress(const FIELD_MATS &fields,
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT_VEC &stress)=0; 
      //* Computes free (T>0)/potential(T=0) energy density 
      //* Units: mvv/L^3 (i.e. for units Real: g/(mol ps^2 A^2) )
      virtual void elastic_energy(const FIELD_MATS &fields,
                                  const GRAD_FIELD_MATS &gradFields,
                                  DENS_MAT &energy) const =0;
      virtual void entropic_energy(const FIELD_MATS &fields,
                                  const GRAD_FIELD_MATS &gradFields,
                                  DENS_MAT &energy) const =0;
      //* Returns the material tangent at a given deformation gradient.
      virtual void tangent(const MATRIX &F, MATRIX &C) const =0;
      //* Creates a linearization for a deformation gradient.
      virtual DENS_VEC elasticity_tensor(const MATRIX &F, MATRIX &C, 
        const ElasticityTensorType type=FIRST_ELASTICITY_TENSOR) const =0;
      virtual DENS_VEC elasticity_tensor(const VECTOR &Fv, MATRIX &C, 
        const ElasticityTensorType type=FIRST_ELASTICITY_TENSOR) const;
    protected:
      double parse_modulus(std::vector<std::string>& line,std::string & units) const;
  };

  /******************************************************************
   * LINEAR MODELS
   ******************************************************************/

  /**
   *  @class  StressLinearElastic
   *  @brief  Class for computing stress for a cubic elastic material
   */ 

  class StressLinearElastic : public Stress
  {
    public:
      StressLinearElastic()
        : c11_(0),c12_(0),c13_(0),c22_(0),c23_(0),c33_(0),
          c44_(0),c55_(0),c66_(0){};
      StressLinearElastic(std::fstream &matfile, 
                          std::string units = "none");
      StressLinearElastic(double c11, double c12, double c44)
        : c11_(c11),c12_(c12),c13_(c12),c22_(c11),c23_(c12),c33_(c11),
          c44_(c44),c55_(c44),c66_(c44){};
      void parse(std::fstream &matfile, std::string units);
      void stress(const FIELD_MATS &fields,
                  const GRAD_FIELD_MATS &gradFields,
                  DENS_MAT_VEC &flux);
      virtual void elastic_energy(const FIELD_MATS &fields,
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT &energy) const;
      virtual void entropic_energy(const FIELD_MATS &fields,
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT &energy) const;
      virtual void tangent(const MATRIX &F, MATRIX &C) const {C=C_;}
      virtual void set_tangent();
      virtual void print_tangent();
      //* Creates a linearization for a deformation gradient.
      virtual DENS_VEC elasticity_tensor(const MATRIX &F, MATRIX &C, 
        const ElasticityTensorType type=FIRST_ELASTICITY_TENSOR) const;
    protected:
      double c11_, c12_, c13_, c22_, c23_, c33_, c44_, c55_, c66_;
      DENS_MAT C_; 
  };

  /**
   *  @class  StressLinearThermoElastic
   *  @brief  Class for computing stress for a cubic elastic material
   */ 

  class StressLinearThermoElastic : public StressLinearElastic
  {
    public:
      StressLinearThermoElastic():alpha1_(0),alpha2_(0),alpha3_(0){};
      StressLinearThermoElastic(std::fstream &matfile, 
                                std::string units = "none");
      StressLinearThermoElastic(double c11, double c12, double c44, 
        double alpha = 0)
        : StressLinearElastic(c11,c12,c44), 
          alpha1_(alpha),alpha2_(alpha),alpha3_(alpha) { }
      void stress(const FIELD_MATS &fields,
                  const GRAD_FIELD_MATS &gradFields,
                  DENS_MAT_VEC &flux);
      virtual void elastic_energy(const FIELD_MATS &fields,
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT &energy) const;
    protected:
      double alpha1_,alpha2_,alpha3_;
      double temperature_;
      DENS_MAT M_;
  };

  /**
   *  @class  StressLinearElasticDamped
   *  @brief  Class for computing stress for a cubic elastic material w/ damping
   */ 

  class StressLinearElasticDamped : public StressLinearElastic
  {
    public:
      StressLinearElasticDamped(std::fstream &matfile, 
                                std::string units = "none");
      StressLinearElasticDamped(double c11, double c12, double c44, 
        double gamma = 0)
        : StressLinearElastic(c11,c12,c44), gamma_(gamma) { }
      void stress(const FIELD_MATS &fields,
                  const GRAD_FIELD_MATS &gradFields,
                  DENS_MAT_VEC &flux);
    protected:
      double gamma_;
  };

  /******************************************************************
   * NON-LINEAR MODELS
   ******************************************************************/

  // forward declarations needed by StressCauchyBorn
  class CbPotential;
  class CBLattice;

  /**
   * Structure of lattice properties needed by StressCauchyBorn.
   */
  struct CbData {
    double e2mvv;           //*> Energy conversion factor (1/mvv2e).
    double boltzmann;       //*> Boltzmann constant (in LAMMPS units)
    double hbar;            //*> Planck's constant (in LAMMPS units)
    double inv_atom_volume; //*> Volume of atom.
    double atom_mass;       //*> Mass of an atom.
    DENS_MAT cell_vectors;  //*> Unit vectors for lattice cells.
    DENS_MAT basis_vectors; //*> Positions of atoms within a lattice cell.
  };

  /**
   * @class StressCauchyBorn
   * @brief Class for computing the stress and elastic constants for a 
   * @brief Cauchy-Born material.
   */


  class StressCauchyBorn : public Stress
  {
    public:
      StressCauchyBorn(std::fstream &matfile, CbData &cb);
      virtual ~StressCauchyBorn(); 
      virtual void initialize(void);
      //* Returns the stress computed from a 0K Cauchy-Born approxmation.
      virtual void stress(const FIELD_MATS &fields, 
                  const GRAD_FIELD_MATS &gradFields,
                  DENS_MAT_VEC &flux);
      //* Computes free (T>0)/potential(T=0) energy density 
      //* Units: mvv/L^3 (i.e. for units Real: g/(mol ps^2 A^2) )
      virtual void elastic_energy(const FIELD_MATS &fields, 
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT &energy) const;
      //* Computes entropic energy density 
      virtual void entropic_energy(const FIELD_MATS &fields, 
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT &energy) const;
      //* Returns the material tangent at a given deformation gradient.
      virtual void tangent(const MATRIX &F, MATRIX &C) const;
      //* Creates a linearization for a deformation gradient.
      virtual DENS_VEC elasticity_tensor(const MATRIX &F, MATRIX &C, 
        const ElasticityTensorType type=FIRST_ELASTICITY_TENSOR) const;
      //* bond stiffness consistent with the einstein freq
      double stiffness() const;
    protected:
      void linearize(MATRIX *F=nullptr);
      CBLattice   *cblattice_;       //*> CbLattice -> makes atom clusters.
      CbPotential *potential_;       //*> CbPotential -> interatomic forces.
      bool makeLinear_;
      StressLinearElastic *cubicMat_; //*> Stores optional linear elastic law.
      bool initialized_;
      double fixed_temperature_;     //*> Specifies a uniform temperature.
      CbData cbdata_;                //*> Lattice & atom volume/mass.
      double alpha_;
      DENS_MAT M_;
  };


    // adaptor to NonLinearSolver
    class ElasticTangentOperator : public TangentOperator {
      public:
        ElasticTangentOperator (Stress * stress,
                                  DENS_VEC & targetP) :
          TangentOperator(),
          stress_(stress),
          targetP_(targetP) {};
        void function(const VECTOR & F, DENS_VEC & R)
        {
          DENS_MAT B;
          tangent(F,R,B); 
        }
        void tangent(const VECTOR & F, DENS_VEC & R, MATRIX & B)
        {
          P_ = stress_->elasticity_tensor(F, B);
          R = P_ - targetP_; 
        }
      private:
        Stress * stress_;
        DENS_VEC targetP_, P_;
    };

    // adaptor to NonLinearSolver
    class SecondElasticTangentOperator : public TangentOperator {
      public:
        SecondElasticTangentOperator (Stress * stress,
                                  DENS_VEC & targetS) :
          TangentOperator(),
          stress_(stress),
          targetS_(targetS) {};
        void function(const VECTOR & U, DENS_VEC & r)
        {
          DENS_MAT C;
          tangent(U,r,C); 
        }
        void tangent(const VECTOR & U, DENS_VEC & r, MATRIX & C)
        {
          S_ = stress_->elasticity_tensor(U, C, SECOND_ELASTICITY_TENSOR);
          
          r = S_ - targetS_; 
        }
      private:
        Stress * stress_;
        DENS_VEC targetS_, S_;
    };
}
#endif 
