#ifndef BOUNDARY_VALUES_HPP
#define BOUNDARY_VALUES_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals.hpp
//  \brief defines BoundaryValues class used for setting BCs on all data types

// C++ headers
#include <string>   // string

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class MeshBlock;
class Hydro;
class Field;
class ParameterInput;
class Radiation;
class Coordinates;
struct FaceField;
struct NeighborBlock;
struct PolarNeighborBlock;

// identifiers for all 6 faces of a MeshBlock
enum BoundaryFace {FACE_UNDEF=-1, INNER_X1=0, OUTER_X1=1, INNER_X2=2, OUTER_X2=3, 
  INNER_X3=4, OUTER_X3=5};

// identifiers for boundary conditions
enum BoundaryFlag {BLOCK_BNDRY=-1, BNDRY_UNDEF=0, REFLECTING_BNDRY=1, OUTFLOW_BNDRY=2,
  USER_BNDRY=3, PERIODIC_BNDRY=4, POLAR_BNDRY=5, POLAR_BNDRY_WEDGE=6, VACUUM_BNDRY=7};

// identifiers for types of neighbor blocks
enum NeighborType {NEIGHBOR_NONE, NEIGHBOR_FACE, NEIGHBOR_EDGE, NEIGHBOR_CORNER};

// identifiers for status of MPI boundary communications
enum BoundaryStatus {BNDRY_WAITING, BNDRY_ARRIVED, BNDRY_COMPLETED};

// flags to mark which variables are reversed across polar boundary
static bool flip_across_pole_hydro[] = {false, false, true, true, false};
static bool flip_across_pole_cr[] = {false, false, true, true};
static bool flip_across_pole_ntc[] = {false, false, true, true};
static bool flip_across_pole_field[] = {false, true, true};

//! \struct NeighborType
//  \brief data to describe MeshBlock neighbors
typedef struct NeighborIndexes {
  int ox1, ox2, ox3, fi1, fi2;
  enum NeighborType type;
  NeighborIndexes() {
    ox1=0; ox2=0; ox3=0; fi1=0; fi2=0;
    type=NEIGHBOR_NONE;
  }
} NeighborIndexes;

//---------------------- prototypes for all BC functions ---------------------------------
void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc, 
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc, 
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc, 
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc, 
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc, 
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc, 
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc, 
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, AthenaArray<Real> &u_ntc,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);



// For radiation boundary condition
void RadReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

void RadOutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadOutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadOutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadOutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadOutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadOutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


void RadVacuumInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadVacuumInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadVacuumInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadVacuumOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadVacuumOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void RadVacuumOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


void RotateHPi_InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            int is, int ie, int js, int je, int ks, int ke);

void RotateHPi_OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            int is, int ie, int js, int je, int ks, int ke);

void RotateHPi_InnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            int is, int ie, int js, int je, int ks, int ke);

void RotateHPi_OuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            int is, int ie, int js, int je, int ks, int ke);


void RotatePi_InnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            int is, int ie, int js, int je, int ks, int ke);

void RotatePi_OuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            int is, int ie, int js, int je, int ks, int ke);



// function to return boundary flag given input string
enum BoundaryFlag GetBoundaryFlag(std::string input_string);


//----------------------------------------------------------------------------------------
//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues {
public:
  BoundaryValues(MeshBlock *pmb, ParameterInput *pin);
  ~BoundaryValues();

  static NeighborIndexes ni[56];
  static int bufid[56];

  void Initialize(void);
  void CheckBoundary(void);
  void StartReceivingForInit(bool cons_and_field);
  void StartReceivingAll(void);
  void ClearBoundaryForInit(bool cons_and_field);
  void ClearBoundaryAll(void);
  void ApplyPhysicalBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
       FaceField &bfdst, AthenaArray<Real> &bcdst, AthenaArray<Real> &crdst,
       AthenaArray<Real> &ntcdst, const Real time, const Real dt);
  void ApplyRadPhysicalBoundaries(AthenaArray<Real> &pdst, const Real time, const Real dt);
  void ProlongateBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
       FaceField &bfdst, AthenaArray<Real> &bcdst, AthenaArray<Real> &rdst, 
       AthenaArray<Real> &crdst, AthenaArray<Real> &ntcdst, const Real time, const Real dt);
  int LoadCellCenteredBoundaryBufferSameLevel(AthenaArray<Real> &src,
                int ns, int ne, Real *buf, const NeighborBlock& nb, int phys);
  int LoadCellCenteredBoundaryBufferToCoarser(AthenaArray<Real> &src,
                          int ns, int ne, Real *buf, AthenaArray<Real> &cbuf,
                                           const NeighborBlock& nb, int phys);
  int LoadCellCenteredBoundaryBufferToFiner(AthenaArray<Real> &src,
                int ns, int ne, Real *buf, const NeighborBlock& nb, int phys);
  void SendCellCenteredBoundaryBuffers(AthenaArray<Real> &src,
          AthenaArray<Real> &src_rad, AthenaArray<Real> &src_cr, 
          AthenaArray<Real> &src_ntc,enum CCBoundaryType type);

  void SetCellCenteredBoundarySameLevel(AthenaArray<Real> &dst,
                          int ns, int ne, Real *buf, const NeighborBlock& nb, 
                                  bool *flip, enum CCBoundaryType type, int &p);
  void SetCellCenteredBoundaryFromCoarser(int ns, int ne,
             Real *buf, AthenaArray<Real> &cbuf,  const NeighborBlock& nb, 
                             bool *flip, enum CCBoundaryType type, int &p);
  void SetCellCenteredBoundaryFromFiner(AthenaArray<Real> &dst,
                          int ns, int ne, Real *buf, const NeighborBlock& nb, 
                                bool *flip, enum CCBoundaryType type, int &p);
  bool ReceiveCellCenteredBoundaryBuffers(AthenaArray<Real> &dst,
                        AthenaArray<Real> &dst_rad, AthenaArray<Real> &dst_cr, 
                        AthenaArray<Real> &dst_ntc, enum CCBoundaryType type);
  void ReceiveCellCenteredBoundaryBuffersWithWait(AthenaArray<Real> &dst,
                        AthenaArray<Real> &dst_rad, AthenaArray<Real> &dst_cr,
                        AthenaArray<Real> &dst_ntc, enum CCBoundaryType type);
  void PolarSingleCellCentered(AthenaArray<Real> &dst, int ns, int ne);
  void PolarSingleRad(AthenaArray<Real> &dst_rad);

  int LoadFieldBoundaryBufferSameLevel(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToCoarser(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToFiner(FaceField &src, Real *buf,
                                     const NeighborBlock& nb);
  void SendFieldBoundaryBuffers(FaceField &src);
  void SetFieldBoundarySameLevel(FaceField &dst, Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromCoarser(Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromFiner(FaceField &dst, Real *buf, const NeighborBlock& nb);
  bool ReceiveFieldBoundaryBuffers(FaceField &dst);
  void ReceiveFieldBoundaryBuffersWithWait(FaceField &dst);
  void PolarSingleField(FaceField &dst);

  void SendFluxCorrection(enum FluxCorrectionType type);
  bool ReceiveFluxCorrection(enum FluxCorrectionType type);

  int LoadEMFBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryPolarBuffer(Real *buf, const PolarNeighborBlock &nb);
  void SendEMFCorrection(void);
  void SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryPolar(Real **buf_list, int num_bufs, bool north);
  void ClearCoarseEMFBoundary(void);
  void AverageEMFBoundary(void);
  void PolarSingleEMF(void);
  bool ReceiveEMFCorrection(void);

private:
  MeshBlock *pmy_block_;  // ptr to MeshBlock containing this BVals

  int nface_, nedge_;
  bool edge_flag_[12];
  int nedge_fine_[12];
  bool firsttime_;

  enum BoundaryStatus hydro_flag_[56], field_flag_[56];
  enum BoundaryStatus flcor_flag_[6][2][2];
  enum BoundaryStatus emfcor_flag_[48];
  enum BoundaryStatus *emf_north_flag_;
  enum BoundaryStatus *emf_south_flag_;
  Real *hydro_send_[56],  *hydro_recv_[56];
  Real *field_send_[56],  *field_recv_[56];
  Real *flcor_send_[6],   *flcor_recv_[6][2][2];
  Real *emfcor_send_[48], *emfcor_recv_[48];
  Real **emf_north_send_, **emf_north_recv_;
  Real **emf_south_send_, **emf_south_recv_;
  AthenaArray<Real> sarea_[2];
  AthenaArray<Real> exc_, rad_exc_;
  int num_north_polar_blocks_, num_south_polar_blocks_;

#ifdef MPI_PARALLEL
  MPI_Request req_hydro_send_[56],  req_hydro_recv_[56];
  MPI_Request req_field_send_[56],  req_field_recv_[56];
  MPI_Request req_flcor_send_[6],   req_flcor_recv_[6][2][2];
  MPI_Request req_emfcor_send_[48], req_emfcor_recv_[48];
  MPI_Request *req_emf_north_send_, *req_emf_north_recv_;
  MPI_Request *req_emf_south_send_, *req_emf_south_recv_;
#endif

  BValFunc_t BoundaryFunction_[6];
  RadBValFunc_t RadBoundaryFunction_[6]; // Function Pointer for radiation

  // temporary
  friend class Mesh;
};

unsigned int CreateBvalsMPITag(int lid, int phys, int bufid);
unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
int BufferID(int dim, bool multilevel);
int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2, int bmax);

#endif // BOUNDARY_VALUES_HPP
