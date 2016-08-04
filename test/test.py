import os
import time
import random
import math
import sys

# 42 is too much of a cliche
random.seed(43)

if len(sys.argv) == 1:
    print "Usage: test.py <executable> [retain]\n"
    print "retain - keep data files\n"
    sys.exit()
lammps = sys.argv[1]

class RandomExperiment():
    def __init__(self, pot):
        self.pot = pot #("BNC", ["B", "N", "C"])
        self.density = 90000. / (60 * 60 * 60 * 8)
        self.ntotal = 600
	self.name = "Random High-Density "+pot[0]
    def get_script(self, pair, pkg, data, suffix):
        data_str = "\n".join(map(lambda d: d.get_commands(suffix), data))
        pot_file, parts = self.pot
	nper = self.ntotal / len(parts)
	create_atoms = "\n".join(map(
	    lambda i: "create_atoms %d random %d %d 1" % (i+1, nper, random.randint(1, 1e9)),
	    range(len(parts))))
	masses = "\n".join(map(
	    lambda i: "mass %d 1.0" % (i+1),
	    range(len(parts))))
        sizes = "-{0:f} {0:f} ".format(math.pow(self.ntotal / self.density / 8, 1./3)) * 3
        return """
        {}
        units metal
	region 1 block {}
        create_box {} 1
        {}
        {}
        pair_style {}
	pair_coeff * * {}.tersoff {}
        fix 1 all nve
        {}
	run 0

        """.format(pkg, sizes, len(parts), create_atoms, masses, pair, pot_file, " ".join(parts), data_str)

class BulkSiStartExperiment():
    def __init__(self):
        self.name = "Start of Bulk Si"
    def get_script(self, pair, pkg, data, suffix):
        data_str = "\n".join(map(lambda d: d.get_commands(suffix), data))
        return """
        {}
        units metal
	atom_style atomic
	lattice diamond 5.431
	region box block 0 20 0 20 0 20
	create_box 1 box
	create_atoms 1 box

	pair_style {}
	mass 1 28.06
	velocity all create 1000.0 {} loop geom
	pair_coeff * * Si.tersoff Si
	neighbor 1.0 bin
	neigh_modify delay 5 every 1
        fix 1 all nve
	timestep 0.001
        {}
	run 0

        """.format(pkg, pair, random.randint(1, 1e9), data_str)


class ForceData():
    def __init__(self, checkers):
        self.checkers = checkers
	self.name = "Force"
    def get_commands(self, suffix):
        return """
	dump dump_force all custom 1 dump{} id x y z fx fy fz
	""".format(suffix)
    def extract_data(self, suffix):
	try:
            with open("dump" + suffix, "r") as f:
                for line in f:
                    for word in line.split():
                        try:
                            yield float(word)
                        except ValueError:
                            pass
	except IOError:
	    print "Warning: Could not open dump file (dump-force-%s)" % suffix
	    pass
   
class GlobalData():
    def __init__(self, var, checkers):
        self.var = var
	self.checkers = checkers
	self.name = "Global " + var
    def get_commands(self, suffix):
        return """
        variable {v}_var equal {v}
	fix {v}_fix all ave/time 1 1 1 v_{v}_var file ave-{v}{s}
	""".format(v = self.var, s = suffix)
    def extract_data(self, suffix):
	try:
            with open("ave-" + self.var + suffix, "r") as f:
                for line in f:
                    for word in line.split():
                        try:
                            yield float(word)
                        except ValueError:
                            pass
	except IOError:
	    print "Warning: Could not open ave file (ave-%s%s)" % (self.var, suffix)
	    pass
   
class Tersoff():
    def __init__(self):
        self.pair = "tersoff"
        self.pair_file = "ters_van"
        self.pkg = ""
	self.name = "tersoff"


class TersoffInlined():
    def __init__(self):
        self.pair = "tersoff/inlined"
        self.pair_file = "ters_inl"
        self.pkg = ""
	self.name = "tersoff/inlined"

class TersoffIntel():
    def __init__(self, run_on, mode):
        self.pair = "tersoff/intel"
        self.pair_file = "ters_int"
        if run_on == "xeon":
            balance = 0
        elif run_on == "phi":
            balance = 1
        else:
            import sys
            print "Error: Unknown device for tersoff/intel style."
            sys.exit(-1)
        self.pkg = """
        package intel 1 mode %s balance %d tptask 1 tpc 1
        package omp 0
        processors * * * grid numa
        """ % (mode, balance)
	self.name = "tersoff/intel %s %s" % (run_on, mode)

class Test():
    def __init__(self):
        import random
        self.trials = []
        self.seed = random.getstate()
        self.cleanup = True
        self.data = [
	    ForceData([MarginOfErrorChecker(0.05)]),
	    GlobalData("evdwl", [MarginOfErrorChecker(0.001)]),
	    GlobalData("press", [MarginOfErrorChecker(0.001)]),
	    ]
    def set_experiment(self, experiment):
        self.experiment = experiment
    def set_reference(self, reference):
        self.reference = reference
    def add_trial(self, trial):
        self.trials.append(trial)
    def set_seed(self, seed):
        self.seed = seed
    def _collect_data(self, of, timestamp):
        import random
        random.setstate(self.seed)
	suffix = "-%15f-%s-tmp" % (timestamp, of.pair_file)
        script_tmp = "script" + suffix
        output_tmp = "output" + suffix
        with open(script_tmp, "w") as f:
            f.write(self.experiment.get_script(of.pair, of.pkg, self.data, suffix))
        os.system('"%s" -in %s > %s' % (lammps, script_tmp, output_tmp))
        return suffix
    def execute(self):
        import time
        run_id = time.time()
        suffix_ref = self._collect_data(self.reference, run_id)
        for trial in self.trials:
	    print "TRIAL    %s vs. %s (%s)" % (self.reference.name, trial.name, self.experiment.name)
            suffix_tr = self._collect_data(trial, run_id)
            for single_data in self.data:
		for checker in single_data.checkers:
	            gen_ref = single_data.extract_data(suffix_ref)
		    gen_tr = single_data.extract_data(suffix_tr)
	            checker.check(gen_ref, gen_tr)
		    checker.report(single_data.name)
            if self.cleanup: os.system("rm *%s" % suffix_tr)
        if self.cleanup: os.system("rm *%s" % suffix_ref)

class FloatChecker():
    def __init__(self):
        pass
    def check(self, gen_ref, gen_trial):
        max_relative_error = 0.
        max_absolute_error = 0.
        found_nan = False
        is_exact = True
        early_end = False
        for elem_ref in gen_ref:
            try:
                elem_trial = gen_trial.next()
                if elem_ref != elem_trial:
                    is_exact = False
                if math.isnan(elem_trial):
                    found_nan = True
                if abs(elem_ref) < 1.e-9:
		    relative_error = abs(elem_trial - elem_ref)
		else:
		    relative_error = abs((elem_ref - elem_trial)/elem_ref)
		absolute_error = abs(elem_ref - elem_trial)
		max_relative_error = max(max_relative_error, relative_error)
		max_absolute_error = max(max_absolute_error, absolute_error)
            except StopIteration:
                early_end = True
                break
        self.max_relative_error = max_relative_error
	self.max_absolute_error = max_absolute_error
	self.found_nan = found_nan
	self.is_exact = is_exact
	self.early_end = early_end
    def report(self, what):
        msg = "%15s. %40s     [%1s%1s%1s %10e,%10e]\033[0m" % (what, self.get_message(), 'E' if self.is_exact else 'I', 'N' if self.found_nan else 'R', 'A' if self.early_end else 'D', self.max_relative_error, self.max_absolute_error)
        if self.is_correct():
            print "\033[32mSUCCESS: " + msg
	else:
            print "\033[91mFAILED : " + msg

class MarginOfErrorChecker(FloatChecker):
    def __init__(self, margin):
        self.margin = margin
    def is_correct(self): return self.max_relative_error < self.margin and not self.early_end and not self.found_nan
    def get_message(self): return "Max. rel. Error (%e)" % self.max_relative_error

class ExactChecker(FloatChecker):
    def is_correct(self): return self.is_exact and not self.early_end and not self.found_nan
    def get_message(self): return "Exact"

if __name__ == "__main__":
    for experiment in [BulkSiStartExperiment(), RandomExperiment(("BNC", ["B", "N", "C"])), RandomExperiment(("GaN", ["Ga", "N"]))]:
        test = Test()
        if len(sys.argv) == 3 and sys.argv[2] == 'retain': test.cleanup = False
        test.set_reference(Tersoff())
        test.experiment = experiment
	for run_on in ["xeon", "phi"]:
	    for mode in ["single", "mixed", "double"]:
                test.add_trial(TersoffIntel(run_on, mode))
        test.execute()
