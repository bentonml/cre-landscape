#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN      # email when started
#SBATCH --mail-user=marylauren.benton@vanderbilt.edu

# Once this job launches, check the node it landed on with:
#      squeue --user=userid
# Then open tunnel from local machine with something like:
#      ssh -L 8888:vmpXXX:7777 vunetid@login.accre.vanderbilt.edu
# This binds port 8888 on localhost to port 7777 on vmpXXX
# Finally, point your local browser to localhost:8888

PORT=7777
echo Execute in a new terminal window:
echo "ssh -L 8888:$HOSTNAME:$PORT $(whoami)@login.accre.vanderbilt.edu"
printf '=%.0s' {1..80}
echo

module load Anaconda3         # load python
source activate enh_gain-loss # load virtual environment

unset XDG_RUNTIME_DIR

# important for removal of .tmp files even after fail/cancel
source /accre/usr/bin/setup_accre_runtime_dir

jt -t solarizedd -f anonymous -fs 12 -nfs 14 -ofs 12 -cellw 90%  # optionally set theme

jupyter notebook --no-browser --ip='*' --port=${PORT}
