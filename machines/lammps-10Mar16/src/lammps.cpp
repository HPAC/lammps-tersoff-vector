/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <string.h>
#include <ctype.h>
#include "lammps.h"
#include "style_angle.h"
#include "style_atom.h"
#include "style_bond.h"
#include "style_command.h"
#include "style_compute.h"
#include "style_dihedral.h"
#include "style_dump.h"
#include "style_fix.h"
#include "style_improper.h"
#include "style_integrate.h"
#include "style_kspace.h"
#include "style_minimize.h"
#include "style_pair.h"
#include "style_region.h"
#include "universe.h"
#include "input.h"
#include "atom.h"
#include "update.h"
#include "neighbor.h"
#include "comm.h"
#include "comm_brick.h"
#include "domain.h"
#include "force.h"
#include "modify.h"
#include "group.h"
#include "output.h"
#include "citeme.h"
#include "accelerator_cuda.h"
#include "accelerator_kokkos.h"
#include "accelerator_omp.h"
#include "accelerator_intel.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   start up LAMMPS
   allocate fundamental classes (memory, error, universe, input)
   parse input switches
   initialize communicators, screen & logfile output
   input is allocated at end after MPI info is setup
------------------------------------------------------------------------- */

LAMMPS::LAMMPS(int narg, char **arg, MPI_Comm communicator)
{
  memory = new Memory(this);
  error = new Error(this);
  universe = new Universe(this,communicator);
  output = NULL;

  screen = NULL;
  logfile = NULL;
  infile = NULL;

  initclock = MPI_Wtime();

  // parse input switches

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int partscreenflag = 0;
  int partlogflag = 0;
  int cudaflag = 0;
  int kokkosflag = 0;
  int restartflag = 0;
  int restartremapflag = 0;
  int citeflag = 1;
  int helpflag = 0;

  suffix = suffix2 = NULL;
  suffix_enable = 0;
  packargs = NULL;
  num_package = 0;
  char *rfile = NULL;
  char *dfile = NULL;
  int wdfirst,wdlast;
  int kkfirst,kklast;

  int npack = 0;
  int *pfirst = NULL;
  int *plast = NULL;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"-partition") == 0 ||
        strcmp(arg[iarg],"-p") == 0) {
      universe->existflag = 1;
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg++;
      while (iarg < narg && arg[iarg][0] != '-') {
        universe->add_world(arg[iarg]);
        iarg++;
      }
    } else if (strcmp(arg[iarg],"-in") == 0 ||
               strcmp(arg[iarg],"-i") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      inflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-screen") == 0 ||
               strcmp(arg[iarg],"-sc") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-log") == 0 ||
               strcmp(arg[iarg],"-l") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-var") == 0 ||
               strcmp(arg[iarg],"-v") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 3;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
    } else if (strcmp(arg[iarg],"-echo") == 0 ||
               strcmp(arg[iarg],"-e") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-pscreen") == 0 ||
               strcmp(arg[iarg],"-ps") == 0) {
      if (iarg+2 > narg)
       error->universe_all(FLERR,"Invalid command-line argument");
      partscreenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-plog") == 0 ||
               strcmp(arg[iarg],"-pl") == 0) {
      if (iarg+2 > narg)
       error->universe_all(FLERR,"Invalid command-line argument");
      partlogflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-cuda") == 0 ||
               strcmp(arg[iarg],"-c") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (strcmp(arg[iarg+1],"on") == 0) cudaflag = 1;
      else if (strcmp(arg[iarg+1],"off") == 0) cudaflag = 0;
      else error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-kokkos") == 0 ||
               strcmp(arg[iarg],"-k") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (strcmp(arg[iarg+1],"on") == 0) kokkosflag = 1;
      else if (strcmp(arg[iarg+1],"off") == 0) kokkosflag = 0;
      else error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
      // delimit any extra args for the Kokkos instantiation
      kkfirst = iarg;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
      kklast = iarg;
    } else if (strcmp(arg[iarg],"-package") == 0 ||
               strcmp(arg[iarg],"-pk") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      memory->grow(pfirst,npack+1,"lammps:pfirst");
      memory->grow(plast,npack+1,"lammps:plast");
      // delimit args for package command invocation
      // any package arg with leading "-" will be followed by numeric digit
      iarg++;
      pfirst[npack] = iarg;
      while (iarg < narg) {
        if (arg[iarg][0] != '-') iarg++;
        else if (isdigit(arg[iarg][1])) iarg++;
        else break;
      }
      plast[npack++] = iarg;
    } else if (strcmp(arg[iarg],"-suffix") == 0 ||
               strcmp(arg[iarg],"-sf") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      delete [] suffix;
      delete [] suffix2;
      suffix2 = NULL;
      suffix_enable = 1;
      // hybrid option to set fall-back for suffix2
      if (strcmp(arg[iarg+1],"hybrid") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
	int n = strlen(arg[iarg+2]) + 1;
	suffix = new char[n];
	strcpy(suffix,arg[iarg+2]);
	n = strlen(arg[iarg+3]) + 1;
	suffix2 = new char[n];
	strcpy(suffix2,arg[iarg+3]);
	iarg += 4;
      } else {
	int n = strlen(arg[iarg+1]) + 1;
	suffix = new char[n];
	strcpy(suffix,arg[iarg+1]);
	iarg += 2;
      }
    } else if (strcmp(arg[iarg],"-reorder") == 0 ||
               strcmp(arg[iarg],"-ro") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (universe->existflag)
        error->universe_all(FLERR,"Cannot use -reorder after -partition");
      universe->reorder(arg[iarg+1],arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"-restart") == 0 ||
               strcmp(arg[iarg],"-r") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      restartflag = 1;
      rfile = arg[iarg+1];
      dfile = arg[iarg+2];
      // check for restart remap flag
      if (strcmp(dfile,"remap") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
        restartremapflag = 1;
        dfile = arg[iarg+3];
        iarg++;
      }
      iarg += 3;
      // delimit any extra args for the write_data command
      wdfirst = iarg;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
      wdlast = iarg;
    } else if (strcmp(arg[iarg],"-nocite") == 0 ||
               strcmp(arg[iarg],"-nc") == 0) {
      citeflag = 0;
      iarg++;
    } else if (strcmp(arg[iarg],"-help") == 0 ||
               strcmp(arg[iarg],"-h") == 0) {
      if (iarg+1 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      helpflag = 1;
      citeflag = 0;
      iarg += 1;
    } else error->universe_all(FLERR,"Invalid command-line argument");
  }

  // if no partition command-line switch, universe is one world with all procs

  if (universe->existflag == 0) universe->add_world(NULL);

  // sum of procs in all worlds must equal total # of procs

  if (!universe->consistent())
    error->universe_all(FLERR,"Processor partitions do not match "
                        "number of allocated processors");

  // universe cannot use stdin for input file

  if (universe->existflag && inflag == 0)
    error->universe_all(FLERR,"Must use -in switch with multiple partitions");

  // if no partition command-line switch, cannot use -pscreen option

  if (universe->existflag == 0 && partscreenflag)
    error->universe_all(FLERR,"Can only use -pscreen with multiple partitions");

  // if no partition command-line switch, cannot use -plog option

  if (universe->existflag == 0 && partlogflag)
    error->universe_all(FLERR,"Can only use -plog with multiple partitions");

  // set universe screen and logfile

  if (universe->me == 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else if (strcmp(arg[screenflag],"none") == 0)
      universe->uscreen = NULL;
    else {
      universe->uscreen = fopen(arg[screenflag],"w");
      if (universe->uscreen == NULL)
        error->universe_one(FLERR,"Cannot open universe screen file");
    }
    if (logflag == 0) {
      if (helpflag == 0) {
        universe->ulogfile = fopen("log.lammps","w");
        if (universe->ulogfile == NULL)
          error->universe_warn(FLERR,"Cannot open log.lammps for writing");
      }
    } else if (strcmp(arg[logflag],"none") == 0)
      universe->ulogfile = NULL;
    else {
      universe->ulogfile = fopen(arg[logflag],"w");
      if (universe->ulogfile == NULL)
        error->universe_one(FLERR,"Cannot open universe log file");
    }
  }

  if (universe->me > 0) {
    if (screenflag == 0) universe->uscreen = stdout;
    else universe->uscreen = NULL;
    universe->ulogfile = NULL;
  }

  // make universe and single world the same, since no partition switch
  // world inherits settings from universe
  // set world screen, logfile, communicator, infile
  // open input script if from file

  if (universe->existflag == 0) {
    screen = universe->uscreen;
    logfile = universe->ulogfile;
    world = universe->uworld;

    if (universe->me == 0) {
      if (inflag == 0) infile = stdin;
      else infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
        char str[128];
        sprintf(str,"Cannot open input script %s",arg[inflag]);
        error->one(FLERR,str);
      }
    }

    if (universe->me == 0) {
      if (screen) fprintf(screen,"LAMMPS (%s)\n",universe->version);
      if (logfile) fprintf(logfile,"LAMMPS (%s)\n",universe->version);
    }

  // universe is one or more worlds, as setup by partition switch
  // split universe communicator into separate world communicators
  // set world screen, logfile, communicator, infile
  // open input script

  } else {
    int me;
    MPI_Comm_split(universe->uworld,universe->iworld,0,&world);
    MPI_Comm_rank(world,&me);

    if (me == 0)
      if (partscreenflag == 0)
       if (screenflag == 0) {
         char str[32];
         sprintf(str,"screen.%d",universe->iworld);
         screen = fopen(str,"w");
         if (screen == NULL) error->one(FLERR,"Cannot open screen file");
       } else if (strcmp(arg[screenflag],"none") == 0)
         screen = NULL;
       else {
         char str[128];
         sprintf(str,"%s.%d",arg[screenflag],universe->iworld);
         screen = fopen(str,"w");
         if (screen == NULL) error->one(FLERR,"Cannot open screen file");
       }
      else if (strcmp(arg[partscreenflag],"none") == 0)
        screen = NULL;
      else {
        char str[128];
        sprintf(str,"%s.%d",arg[partscreenflag],universe->iworld);
        screen = fopen(str,"w");
        if (screen == NULL) error->one(FLERR,"Cannot open screen file");
      } else screen = NULL;

    if (me == 0)
      if (partlogflag == 0)
       if (logflag == 0) {
         char str[32];
         sprintf(str,"log.lammps.%d",universe->iworld);
         logfile = fopen(str,"w");
         if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
       } else if (strcmp(arg[logflag],"none") == 0)
         logfile = NULL;
       else {
         char str[128];
         sprintf(str,"%s.%d",arg[logflag],universe->iworld);
         logfile = fopen(str,"w");
         if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
       }
      else if (strcmp(arg[partlogflag],"none") == 0)
        logfile = NULL;
      else {
        char str[128];
        sprintf(str,"%s.%d",arg[partlogflag],universe->iworld);
        logfile = fopen(str,"w");
        if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
      } else logfile = NULL;

    if (me == 0) {
      infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
        char str[128];
        sprintf(str,"Cannot open input script %s",arg[inflag]);
        error->one(FLERR,str);
      }
    } else infile = NULL;

    // screen and logfile messages for universe and world

    if (universe->me == 0) {
      if (universe->uscreen) {
        fprintf(universe->uscreen,"LAMMPS (%s)\n",universe->version);
        fprintf(universe->uscreen,"Running on %d partitions of processors\n",
                universe->nworlds);
      }
      if (universe->ulogfile) {
        fprintf(universe->ulogfile,"LAMMPS (%s)\n",universe->version);
        fprintf(universe->ulogfile,"Running on %d partitions of processors\n",
                universe->nworlds);
      }
    }

    if (me == 0) {
      if (screen) {
        fprintf(screen,"LAMMPS (%s)\n",universe->version);
        fprintf(screen,"Processor partition = %d\n",universe->iworld);
      }
      if (logfile) {
        fprintf(logfile,"LAMMPS (%s)\n",universe->version);
        fprintf(logfile,"Processor partition = %d\n",universe->iworld);
      }
    }
  }

  // check consistency of datatype settings in lmptype.h

  if (sizeof(smallint) != sizeof(int))
    error->all(FLERR,"Smallint setting in lmptype.h is invalid");
  if (sizeof(imageint) < sizeof(smallint))
    error->all(FLERR,"Imageint setting in lmptype.h is invalid");
  if (sizeof(tagint) < sizeof(smallint))
    error->all(FLERR,"Tagint setting in lmptype.h is invalid");
  if (sizeof(bigint) < sizeof(imageint) || sizeof(bigint) < sizeof(tagint))
    error->all(FLERR,"Bigint setting in lmptype.h is invalid");

  int mpisize;
  MPI_Type_size(MPI_LMP_TAGINT,&mpisize);
  if (mpisize != sizeof(tagint))
      error->all(FLERR,"MPI_LMP_TAGINT and tagint in "
                 "lmptype.h are not compatible");
  MPI_Type_size(MPI_LMP_BIGINT,&mpisize);
  if (mpisize != sizeof(bigint))
      error->all(FLERR,"MPI_LMP_BIGINT and bigint in "
                 "lmptype.h are not compatible");

#ifdef LAMMPS_SMALLBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 ||
      sizeof(tagint) != 4 || sizeof(bigint) != 8)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif
#ifdef LAMMPS_BIGBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 8 ||
      sizeof(tagint) != 8 || sizeof(bigint) != 8)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif
#ifdef LAMMPS_SMALLSMALL
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 ||
      sizeof(tagint) != 4 || sizeof(bigint) != 4)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif

  // error check on accelerator packages

  if (cudaflag == 1 && kokkosflag == 1)
    error->all(FLERR,"Cannot use -cuda on and -kokkos on together");

  // create Cuda class if USER-CUDA installed, unless explicitly switched off
  // instantiation creates dummy Cuda class if USER-CUDA is not installed

  cuda = NULL;
  if (cudaflag == 1) {
    cuda = new Cuda(this);
    if (!cuda->cuda_exists)
      error->all(FLERR,"Cannot use -cuda on without USER-CUDA installed");
  }

  int me;
  MPI_Comm_rank(world,&me);
  if (cuda && me == 0) error->message(FLERR,"USER-CUDA mode is enabled");

  // create Kokkos class if KOKKOS installed, unless explicitly switched off
  // instantiation creates dummy Kokkos class if KOKKOS is not installed
  // add args between kkfirst and kklast to Kokkos instantiation

  kokkos = NULL;
  if (kokkosflag == 1) {
    kokkos = new KokkosLMP(this,kklast-kkfirst,&arg[kkfirst]);
    if (!kokkos->kokkos_exists)
      error->all(FLERR,"Cannot use -kokkos on without KOKKOS installed");
  }

  // allocate CiteMe class if enabled

  if (citeflag) citeme = new CiteMe(this);
  else citeme = NULL;

  // allocate input class now that MPI is fully setup

  input = new Input(this,narg,arg);

  // copy package cmdline arguments
  if (npack > 0) {
    num_package = npack;
    packargs = new char**[npack];
    for (int i=0; i < npack; ++i) {
      int n = plast[i] - pfirst[i];
      packargs[i] = new char*[n+1];
      for (int j=0; j < n; ++j)
        packargs[i][j] = strdup(arg[pfirst[i]+j]);
      packargs[i][n] = NULL;
    }
    memory->destroy(pfirst);
    memory->destroy(plast);
  }

  // allocate top-level classes

  create();
  post_create();

  // if helpflag set, print help and quit with "success" status

  if (helpflag) {
    if (universe->me == 0 && screen) help();
    error->done(0);
  }

  // if restartflag set, invoke 2 commands and quit
  // add args between wdfirst and wdlast to write_data command
  // also add "noinit" to prevent write_data from doing system init

  if (restartflag) {
    char cmd[128];
    sprintf(cmd,"read_restart %s\n",rfile);
    if (restartremapflag) strcat(cmd," remap\n");
    input->one(cmd);
    sprintf(cmd,"write_data %s",dfile);
    for (iarg = wdfirst; iarg < wdlast; iarg++)
      sprintf(&cmd[strlen(cmd)]," %s",arg[iarg]);
    strcat(cmd," noinit\n");
    input->one(cmd);
    error->done(0);
  }
}

/* ----------------------------------------------------------------------
   shutdown LAMMPS
   delete top-level classes
   close screen and log files in world and universe
   output files were already closed in destroy()
   delete fundamental classes
------------------------------------------------------------------------- */

LAMMPS::~LAMMPS()
{
  const int me = comm->me;

  destroy();
  delete citeme;

  if (num_package) {
    for (int i = 0; i < num_package; i++) {
      for (char **ptr = packargs[i]; *ptr != NULL; ++ptr)
        free(*ptr);
      delete[] packargs[i];
    }
    delete[] packargs;
  }
  num_package = 0;
  packargs = NULL;

  double totalclock = MPI_Wtime() - initclock;
  if ((me == 0) && (screen || logfile)) {
    char outtime[128];
    int seconds = fmod(totalclock,60.0);
    totalclock  = (totalclock - seconds) / 60.0;
    int minutes = fmod(totalclock,60.0);
    int hours = (totalclock - minutes) / 60.0;
    sprintf(outtime,"Total wall time: "
            "%d:%02d:%02d\n", hours, minutes, seconds);
    if (screen) fputs(outtime,screen);
    if (logfile) fputs(outtime,logfile);
  }

  if (universe->nworlds == 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    logfile = NULL;
    if (screen != stdout) screen = NULL;
  } else {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    if (universe->ulogfile) fclose(universe->ulogfile);
    logfile = NULL;
    if (screen != stdout) screen = NULL;
  }

  if (infile && infile != stdin) fclose(infile);

  if (world != universe->uworld) MPI_Comm_free(&world);

  delete cuda;
  delete kokkos;
  delete [] suffix;
  delete [] suffix2;

  delete input;
  delete universe;
  delete error;
  delete memory;
}

/* ----------------------------------------------------------------------
   allocate single instance of top-level classes
   fundamental classes are allocated in constructor
   some classes have package variants
------------------------------------------------------------------------- */

void LAMMPS::create()
{
  // Comm class must be created before Atom class
  // so that nthreads is defined when create_avec invokes grow()

  if (cuda) comm = new CommCuda(this);
  else if (kokkos) comm = new CommKokkos(this);
  else comm = new CommBrick(this);

  if (cuda) neighbor = new NeighborCuda(this);
  else if (kokkos) neighbor = new NeighborKokkos(this);
  else neighbor = new Neighbor(this);

  if (cuda) domain = new DomainCuda(this);
  else if (kokkos) domain = new DomainKokkos(this);
#ifdef LMP_USER_OMP
  else domain = new DomainOMP(this);
#else
  else domain = new Domain(this);
#endif

  if (kokkos) atom = new AtomKokkos(this);
  else atom = new Atom(this);
  atom->create_avec("atomic",0,NULL,1);

  group = new Group(this);
  force = new Force(this);    // must be after group, to create temperature

  if (cuda) modify = new ModifyCuda(this);
  else if (kokkos) modify = new ModifyKokkos(this);
  else modify = new Modify(this);

  output = new Output(this);  // must be after group, so "all" exists
                              // must be after modify so can create Computes
  update = new Update(this);  // must be after output, force, neighbor
  timer = new Timer(this);
}

/* ----------------------------------------------------------------------
   check suffix consistency with installed packages
   invoke package-specific deafult package commands
     only invoke if suffix is set and enabled
     also check if suffix2 is set
   called from LAMMPS constructor and after clear() command
     so that package-specific core classes have been instantiated
------------------------------------------------------------------------- */

void LAMMPS::post_create()
{
  // default package commands triggered by "-c on" and "-k on"

  if (cuda && cuda->cuda_exists) input->one("package cuda 1");
  if (kokkos && kokkos->kokkos_exists) input->one("package kokkos");

  // suffix will always be set if suffix_enable = 1
  // check that USER-CUDA and KOKKOS package classes were instantiated
  // check that GPU, INTEL, USER-OMP fixes were compiled with LAMMPS

  if (!suffix_enable) return;

  if (strcmp(suffix,"cuda") == 0 && (cuda == NULL || cuda->cuda_exists == 0))
    error->all(FLERR,"Using suffix cuda without USER-CUDA package enabled");
  if (strcmp(suffix,"gpu") == 0 && !modify->check_package("GPU"))
    error->all(FLERR,"Using suffix gpu without GPU package installed");
  if (strcmp(suffix,"intel") == 0 && !modify->check_package("INTEL"))
    error->all(FLERR,"Using suffix intel without USER-INTEL package installed");
  if (strcmp(suffix,"kk") == 0 &&
      (kokkos == NULL || kokkos->kokkos_exists == 0))
    error->all(FLERR,"Using suffix kk without KOKKOS package enabled");
  if (strcmp(suffix,"omp") == 0 && !modify->check_package("OMP"))
    error->all(FLERR,"Using suffix omp without USER-OMP package installed");

  if (strcmp(suffix,"gpu") == 0) input->one("package gpu 1");
  if (strcmp(suffix,"intel") == 0) input->one("package intel 1");
  if (strcmp(suffix,"omp") == 0) input->one("package omp 0");

  if (suffix2) {
    if (strcmp(suffix2,"gpu") == 0) input->one("package gpu 1");
    if (strcmp(suffix2,"intel") == 0) input->one("package intel 1");
    if (strcmp(suffix2,"omp") == 0) input->one("package omp 0");
  }

  // invoke any command-line package commands

  if (num_package) {
    char str[256];
    for (int i = 0; i < num_package; i++) {
      strcpy(str,"package");
      for (char **ptr = packargs[i]; *ptr != NULL; ++ptr) {
        if (strlen(str) + strlen(*ptr) + 2 > 256)
          error->all(FLERR,"Too many -pk arguments in command line");
        strcat(str," ");
        strcat(str,*ptr);
      }
      input->one(str);
    }
  }
}

/* ----------------------------------------------------------------------
   initialize top-level classes
   do not initialize Timer class, other classes like Run() do that explicitly
------------------------------------------------------------------------- */

void LAMMPS::init()
{
  update->init();
  force->init();         // pair must come after update due to minimizer
  domain->init();
  atom->init();          // atom must come after force and domain
                         //   atom deletes extra array
                         //   used by fix shear_history::unpack_restart()
                         //   when force->pair->gran_history creates fix ??
                         //   atom_vec init uses deform_vremap
  modify->init();        // modify must come after update, force, atom, domain
  neighbor->init();      // neighbor must come after force, modify
  comm->init();          // comm must come after force, modify, neighbor, atom
  output->init();        // output must come after domain, force, modify
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
------------------------------------------------------------------------- */

void LAMMPS::destroy()
{
  delete update;
  update = NULL;

  delete neighbor;
  neighbor = NULL;

  delete comm;
  comm = NULL;

  delete force;
  force = NULL;

  delete group;
  group = NULL;

  delete output;
  output = NULL;

  delete modify;          // modify must come after output, force, update
                          //   since they delete fixes
  modify = NULL;

  delete domain;          // domain must come after modify
                          //   since fix destructors access domain
  domain = NULL;

  delete atom;            // atom must come after modify, neighbor
                          //   since fixes delete callbacks in atom
  atom = NULL;

  delete timer;
  timer = NULL;
}

/* ----------------------------------------------------------------------
   help message for command line options and styles present in executable
------------------------------------------------------------------------- */

void LAMMPS::help()
{
  fprintf(screen,
          "\nCommand line options:\n\n"
          "-cuda on/off                : turn CUDA mode on or off (-c)\n"
          "-echo none/screen/log/both  : echoing of input script (-e)\n"
          "-help                       : print this help message (-h)\n"
          "-in filename                : read input from file, not stdin (-i)\n"
          "-kokkos on/off ...          : turn KOKKOS mode on or off (-k)\n"
          "-log none/filename          : where to send log output (-l)\n"
          "-nocite                     : disable writing log.cite file (-nc)\n"
          "-package style ...          : invoke package command (-pk)\n"
          "-partition size1 size2 ...  : assign partition sizes (-p)\n"
          "-plog basename              : basename for partition logs (-pl)\n"
          "-pscreen basename           : basename for partition screens (-ps)\n"
          "-restart rfile dfile ...    : convert restart to data file (-r)\n"
          "-reorder topology-specs     : processor reordering (-r)\n"
          "-screen none/filename       : where to send screen output (-sc)\n"
          "-suffix cuda/gpu/opt/omp    : style suffix to apply (-sf)\n"
          "-var varname value          : set index style variable (-v)\n\n");

  fprintf(screen,"Style options compiled with this executable\n\n");

  int pos = 80;
  fprintf(screen,"* Atom styles:\n");
#define ATOM_CLASS
#define AtomStyle(key,Class) print_style(#key,pos);
#include "style_atom.h"
#undef ATOM_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Integrate styles:\n");
#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) print_style(#key,pos);
#include "style_integrate.h"
#undef INTEGRATE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Minimize styles:\n");
#define MINIMIZE_CLASS
#define MinimizeStyle(key,Class) print_style(#key,pos);
#include "style_minimize.h"
#undef MINIMIZE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Pair styles:\n");
#define PAIR_CLASS
#define PairStyle(key,Class) print_style(#key,pos);
#include "style_pair.h"
#undef PAIR_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Bond styles:\n");
#define BOND_CLASS
#define BondStyle(key,Class) print_style(#key,pos);
#include "style_bond.h"
#undef BOND_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Angle styles:\n");
#define ANGLE_CLASS
#define AngleStyle(key,Class) print_style(#key,pos);
#include "style_angle.h"
#undef ANGLE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Dihedral styles:\n");
#define DIHEDRAL_CLASS
#define DihedralStyle(key,Class) print_style(#key,pos);
#include "style_dihedral.h"
#undef DIHEDRAL_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Improper styles:\n");
#define IMPROPER_CLASS
#define ImproperStyle(key,Class) print_style(#key,pos);
#include "style_improper.h"
#undef IMPROPER_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* KSpace styles:\n");
#define KSPACE_CLASS
#define KSpaceStyle(key,Class) print_style(#key,pos);
#include "style_kspace.h"
#undef KSPACE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Fix styles\n");
#define FIX_CLASS
#define FixStyle(key,Class) print_style(#key,pos);
#include "style_fix.h"
#undef FIX_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Compute styles:\n");
#define COMPUTE_CLASS
#define ComputeStyle(key,Class) print_style(#key,pos);
#include "style_compute.h"
#undef COMPUTE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Region styles:\n");
#define REGION_CLASS
#define RegionStyle(key,Class) print_style(#key,pos);
#include "style_region.h"
#undef REGION_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Dump styles:\n");
#define DUMP_CLASS
#define DumpStyle(key,Class) print_style(#key,pos);
#include "style_dump.h"
#undef DUMP_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Command styles\n");
#define COMMAND_CLASS
#define CommandStyle(key,Class) print_style(#key,pos);
#include "style_command.h"
#undef COMMAND_CLASS
  fprintf(screen,"\n");
}

/* ----------------------------------------------------------------------
   print style names in columns
   skip any style that starts with upper-case letter, since internal
------------------------------------------------------------------------- */

void LAMMPS::print_style(const char *str, int &pos)
{
  if (isupper(str[0])) return;

  int len = strlen(str);
  if (pos+len > 80) {
    fprintf(screen,"\n");
    pos = 0;
  }

  if (len < 16) {
    fprintf(screen,"%-16s",str);
    pos += 16;
  } else if (len < 32) {
    fprintf(screen,"%-32s",str);
    pos += 32;
  } else if (len < 48) {
    fprintf(screen,"%-48s",str);
    pos += 48;
  } else if (len < 64) {
    fprintf(screen,"%-64s",str);
    pos += 64;
  } else {
    fprintf(screen,"%-80s",str);
    pos += 80;
  }
}
