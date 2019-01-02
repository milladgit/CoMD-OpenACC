
import os, sys

job_template = '''#!/bin/bash

### Set the job name
#SBATCH -J %s
#SBATCH --mail-user=mghane2@uh.edu
#SBATCH --mail-type=begin    # email me when the job starts
#SBATCH --mail-type=end      # email me when the job finishes
#SBATCH -o %s
#SBATCH -e %s

#SBATCH -p hsw_p100

### Specify the number of cpus for your job.
#SBATCH -N 1                 # total number of nodes
#SBATCH --exclusive


source /etc/profile.d/modules.sh


module load cuda pgi

export COMD_GANG_COUNT=%d
export COMD_VECTOR_COUNT=%d

%s 

'''



app = "CoMD-acc-v09-tesla"
cmd = "nvprof --system-profiling on %s " % (app)

count = 0

for gang_size in range(200, 2001, 100):
	for vector_len in range(16, 129, 16):
		name = "%d-%d" % (gang_size, vector_len)

		output_file = "output-%s" % (name)
		error_file = "error-%s" % (name)

		job = job_template % (name, output_file, error_file, gang_size, vector_len, cmd)

		filename = "%s.job" % (name)

		f = open(filename, "w")
		f.write(job)
		f.close()

		os.system("sbatch %s" % (filename))

