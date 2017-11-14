# For easy data analysis, store results in a sqlite database and a dat file,
# such that we can do data analysis with both SQL and R.
# This tool will automatically discover results in the out/ subdirectory.
# It may be necessary to remove data.db/data.dat between subsequent runs.

import sqlite3
import glob
import re

conn = sqlite3.connect('data.db')
ff = open("data.dat", "w")
ff.write("hostname precision arch vect loop run_index sim runtime\n")

c = conn.cursor()

c.execute('CREATE TABLE runs (id INTEGER PRIMARY KEY ASC, hostname VARCHAR(255), sim VARCHAR(255), precision VARCHAR(31), arch CHARACTER(1), vect CHARACTER(1), loop CHARACTER(1), run_index INTEGER, runtime DOUBLE, lastline VARCHAR(255))')

re_filename = re.compile("out/(?P<hostname>.*)\\.(?P<sim>[^.]*)\\.(?P<precision>(single|double|mixed|original))\\.lmp_intel_eval_(?P<arch>.)(?P<vect>(s|v))(?P<loop>(i|j))\\.(?P<run_index>.)")
re_loop = re.compile("Loop time of *(?P<runtime>[0-9]*\\.[0-9]*) on")

for fn in glob.glob("out/*"):
    # parse filename
    m = re_filename.match(fn)
    info = m.groupdict()
    # parse file
    with open(fn, "r") as f:
        for line in f:
            m = re_loop.match(line)
            if m:
                info.update(m.groupdict())
                c.execute('INSERT INTO runs (hostname, precision, arch, vect, loop, run_index, sim, runtime, lastline) VALUES (:hostname, :precision, :arch, :vect, :loop, :run_index, :sim, :runtime, :lastline)', info)
                ff.write("{hostname} {precision} {arch} {vect} {loop} {run_index} {sim} {runtime}\n".format(**info))
                break
            info["lastline"] = line

ff.close()
    
c.execute('select hostname, arch, precision, vect, min(timesum) from (select runid, sum(time) as timesum from runtimes group by runid) INNER JOIN runs on runs.id = runid group by hostname, arch, precision, vect order by hostname, arch, length(precision) desc, precision, vect desc')
lasthostname = ""
for row in c.fetchall():
    hostname = row[0]
    arch = row[1]
    precision = row[2]
    vect = row[3]
    timesum = row[4]
    if lasthostname == hostname:
        hostname = ""
    else:
       lasthostname = hostname
    if vect == "v":
        print "%20s %s %10s %s %f" % (hostname, arch, precision, vect, timesum),
    else:
        print "%f" % timesum
    if precision == "original": print
    #print "%20s %s %10s %s %f" % (hostname, arch, precision, vect, timesum)

conn.commit()
conn.close()
