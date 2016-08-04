// -*- c++ -*-

#ifndef COLVARPROXY_LAMMPS_H
#define COLVARPROXY_LAMMPS_H

#include "colvarmodule.h"
#include "colvarproxy.h"

#include "lammps.h"
#include "domain.h"
#include "force.h"
#include "random_park.h"
#include "update.h"

#include <string>
#include <vector>
#include <iostream>

#ifndef COLVARPROXY_VERSION
#define COLVARPROXY_VERSION "2016-02-28"
#endif

/* struct for packed data communication of coordinates and forces. */
struct commdata {
  int tag,type;
  double x,y,z,m,q;
};

inline std::ostream & operator<< (std::ostream &out, const commdata &cd)
{
  out << " (" << cd.tag << "/" << cd.type << ": "
      << cd.x << ", " << cd.y << ", " << cd.z << ") ";
  return out;
}

/// \brief Communication between colvars and LAMMPS
/// (implementation of \link colvarproxy \endlink)
class colvarproxy_lammps : public colvarproxy {

  // LAMMPS specific data objects and flags
 protected:

  // pointers to LAMMPS class instances
  class LAMMPS_NS::LAMMPS *_lmp;
  class LAMMPS_NS::RanPark *_random;

  // state of LAMMPS properties
  double t_target, my_timestep, my_boltzmann, my_angstrom;
  double bias_energy;
  int  restart_every;
  int  previous_step;

  bool first_timestep;
  bool system_force_requested;
  bool do_exit;

  // std::vector<int>          colvars_atoms;
  // std::vector<size_t>       colvars_atoms_ncopies;
  // std::vector<struct commdata> positions;
  // std::vector<struct commdata> total_forces;
  // std::vector<struct commdata> applied_forces;
  // std::vector<struct commdata> previous_applied_forces;

  std::vector<int>          atoms_types;

  MPI_Comm inter_comm;     // MPI comm with 1 root proc from each world
  int inter_me, inter_num; // rank for the inter replica comm

 public:
  friend class cvm::atom;
  colvarproxy_lammps(LAMMPS_NS::LAMMPS *lmp, const char *,
                     const char *, const int, const double, MPI_Comm);
  virtual ~colvarproxy_lammps();
  void init(const char*);
  int setup();

 // disable default and copy constructor
 private:
  colvarproxy_lammps() {};
  colvarproxy_lammps(const colvarproxy_lammps &) {};

  // methods for lammps to move data or trigger actions in the proxy
 public:
  void set_temperature(double t) { t_target = t; };
  bool need_system_forces() const { return  system_force_requested; };
  bool want_exit() const { return do_exit; };

  // perform colvars computation. returns biasing energy
  double compute();

  // dump status to string
  void serialize_status(std::string &);

  // set status from string
  bool deserialize_status(std::string &);

  // implementation of pure methods from base class
 public:

  inline cvm::real unit_angstrom() { return my_angstrom; };
  inline cvm::real boltzmann() { return my_boltzmann; };
  inline cvm::real temperature() { return t_target; };
  inline cvm::real dt() { return my_timestep; }; // return _lmp->update->dt * _lmp->force->femtosecond; };

  inline size_t restart_frequency() { return restart_every; };

  void add_energy(cvm::real energy) { bias_energy += energy; };
  void request_system_force(bool yesno) { system_force_requested = yesno; };

  void log(std::string const &message);
  void error(std::string const &message);
  void fatal_error(std::string const &message);
  void exit(std::string const &message);

  cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                 cvm::atom_pos const &pos2);
  cvm::real position_dist2(cvm::atom_pos const &pos1,
                           cvm::atom_pos const &pos2);
  void select_closest_image(cvm::atom_pos &pos,
                            cvm::atom_pos const &ref_pos);

  int backup_file(char const *filename);

  cvm::real rand_gaussian(void) { return _random->gaussian(); };

  int init_atom(int atom_number);
  int check_atom_id(int atom_number);

  inline std::vector<int> *modify_atom_types() { return &atoms_types; }

  // implementation of optional methods from base class
 public:
  // Multi-replica support
  // Indicate if multi-replica support is available and active
  virtual bool replica_enabled() { return (inter_comm != MPI_COMM_NULL); }

  // Index of this replica
  virtual int replica_index() { return inter_me; }

  // Total number of replica
  virtual int replica_num() { return inter_num; }

  // Synchronize replica
  virtual void replica_comm_barrier();

  // Receive data from other replica
  virtual int replica_comm_recv(char* msg_data, int buf_len, int src_rep);

  // Send data to other replica
  virtual int replica_comm_send(char* msg_data, int msg_len, int dest_rep);
};

#endif

