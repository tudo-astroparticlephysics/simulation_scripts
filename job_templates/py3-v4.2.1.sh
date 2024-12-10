#!/bin/bash
#PBS -l nodes=1:ppn={cpus}
#PBS -l pmem={memory}
#PBS -l mem={memory}
#PBS -l vmem={memory}
#PBS -l pvmem={memory}
#PBS -l walltime={walltime}
#PBS -o {processing_folder}/logs/{step_name}_run_{run_number}_${PBS_JOBID}.out
#PBS -e {processing_folder}/logs/{step_name}_run_{run_number}_${PBS_JOBID}.err
#PBS -q long
#PBS -S /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/icetray-start
FINAL_OUT={final_out}
KEEP_CRASHED_FILES={keep_crashed_files}


cd
echo 'PWD: '$(pwd)
echo 'Starting job on Host: '$HOSTNAME
echo 'Loading py3-v4.2.1'
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/setup.sh`

# optionally set a different python userbase
# Note: this should be avoided if possible
export PYTHONUSERBASE={python_user_base}
echo 'Using PYTHONUSERBASE: '${PYTHONUSERBASE}

export ENV_SITE_PACKAGES=$(find ${PYTHONUSERBASE}/lib* -maxdepth 2 -type d -name "site-packages")
export PYTHONPATH=$ENV_SITE_PACKAGES:$PYTHONPATH
export PATH=$PYTHONUSERBASE/bin:$PATH
echo 'Using PYTHONPATH: '${PYTHONPATH}

echo $FINAL_OUT
if [ -z ${PBS_JOBID} ] && [ -z ${_CONDOR_SCRATCH_DIR} ]
then
    echo 'Running Script w/o temporary scratch'
    {script_folder}/steps/{step_name}.py {yaml_copy} {run_number} --no-scratch
    ICETRAY_RC=$?
    echo 'IceTray finished with Exit Code: ' $ICETRAY_RC
    if [ $ICETRAY_RC -ne 0 ] && [ $KEEP_CRASHED_FILES -eq 0 ] ; then
        echo 'Deleting partially processed file! ' $FINAL_OUT
        rm $FINAL_OUT
    fi
else
    echo 'Running Script w/ temporary scratch'
    if [ -z ${_CONDOR_SCRATCH_DIR} ]
    then
        cd /scratch/${USER}
    else
        cd ${_CONDOR_SCRATCH_DIR}
    fi
    {script_folder}/steps/{step_name}.py {yaml_copy} {run_number} --scratch
    ICETRAY_RC=$?
    echo 'IceTray finished with Exit Code: ' $ICETRAY_RC
    if [ $ICETRAY_RC -eq 0 ] || [ $KEEP_CRASHED_FILES -eq 1 ]; then
        cp *.i3.bz2 {output_folder}
    fi
    rm *.i3.bz2
fi
exit $ICETRAY_RC

