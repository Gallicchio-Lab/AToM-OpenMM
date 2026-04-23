#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --output=<JOBNAME>.log
#SBATCH --error=<JOBNAME>.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --chdir=<WORKDIR>
#SBATCH -t 26:00:00

# term handler
# the function is executed once the job gets the SIGTERM signal
# because of preemption (or scancel, do scancel -s 9 jobid to
# force hard termination)
term_handler()
{
    #send sigterm to AToM to have it checkpoint
    echo "executing term_handler at $(date)"
    if [ -n "$ATOMPID" ] ; then
        echo "Sending SIGTERM signal to process $ATOMPID "
        kill -SIGTERM $ATOMPID
    fi

    echo "Waiting 200 secs for the signal to take effect"
    sleep 200
}
# declare the function handling the TERM signal
if [ -n "$SLURM_JOB_ID" ]; then
    trap 'term_handler' TERM
fi

jobname=<JOBNAME>
sdir=<SCRIPTS_DIR>

. <CONDADIR>/../../bin/activate <CONDAENV>

echo "Running on $(hostname)"

CMD=(
    python -u ${sdir}/run-atm.py
    --optionsYAMLinFile ${sdir}/defaults.yaml
    --jobBasename ${jobname}
    --receptorinFile ${sdir}/../receptor/<RCPT>.pdb
    --LIG1inFile ${sdir}/../ligands/<LIG1>.sdf
    --plotOutFile ${jobname}.png
)

if [ -n "$SLURM_JOB_ID" ]; then
    "${CMD[@]}" &
    ATOMPID=$!
    wait "$ATOMPID"
else
    "${CMD[@]}"
fi
