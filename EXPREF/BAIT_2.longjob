#!/bin/bash -l  
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job.
#SBATCH --export=ALL
# Define job name
#SBATCH -J BAIT_2
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o nemo.%u.%N.%j.out
# Define a standard error file
#SBATCH -e nemo.%u.%N.%j.err
# Request the partition
#SBATCH -p nodes
# Request the number of nodes
#SBATCH -N 2
# Request the total number of cores
#SBATCH -n 80 
# This asks for 0 days, 8 hour, 0 minutes and 0 seconds of time.
#SBATCH -t 2-23:59:59
# Specify memory per core
#SBATCH --mem-per-cpu=9000M
#
#
# define some information about the job
ulimit -s unlimited
export NPROC=80
export XPROC=0
let RUNPROC=$NPROC-$XPROC
export ENAM=BAIT_2
export YEAR=2200
export CONT=1
export ENDYR=3000
export STYR=${YEAR}
#export RBLD="/users/atagliab/dev_r8832_PISCO/NEMOGCM/TOOLS/REBUILD_NEMO"
#export RBLD2="/mnt/data2/users/atagliab/NEMOv4.0/tools/REBUILD"
export RBLD2="/users/atagliab/RUNS/NEMO42/BAIT/${ENAM}"
#export RBLD2="/users/atagliab/NEMOv4.0/tools/REBUILD"
export BIN="/users/atagliab/RUNS/NEMO42/BAIT/${ENAM}"
export XBIN="/mnt/data2/users/atagliab/xios-2.5i/bin"
export OUTDIR="/mnt/data2/users/atagliab/NEMO_OUT/NEMO42/BAIT/${ENAM}"
export RUNDIR="/mnt/data2/users/atagliab/RUNDIR4"
export RUN="${OUTDIR}/OUT"
#export RUN="/mnt/lustre/users/atagliab/${ENAM}"
echo =========================================================
echo Job submitted date = `date`
date_start=`date +%s`

hostname

echo "Print the following environmetal variables:"
echo "Job name                     : $SLURM_JOB_NAME"
echo "Job ID                       : $SLURM_JOB_ID"
echo "Job user                     : $SLURM_JOB_USER"
echo "Job array index              : $SLURM_ARRAY_TASK_ID"
echo "Submit directory             : $SLURM_SUBMIT_DIR"
echo "Temporary directory          : $TMPDIR"
echo "Submit host                  : $SLURM_SUBMIT_HOST"
echo "Queue/Partition name         : $SLURM_JOB_PARTITION"
echo "Node list                    : $SLURM_JOB_NODELIST"
echo "Hostname of 1st node         : $HOSTNAME"
echo "Number of nodes allocated    : $SLURM_JOB_NUM_NODES or $SLURM_NNODES"
echo "Number of processes          : $SLURM_NTASKS"
echo "Number of processes per node : $SLURM_TASKS_PER_NODE"
echo "Requested tasks per node     : $SLURM_NTASKS_PER_NODE"
echo "Requested CPUs per task      : $SLURM_CPUS_PER_TASK"
echo "Scheduling priority          : $SLURM_PRIO_PROCESS"

cd ${RUN}

echo 'In directory: ' $RUN

if [ $CONT -eq 1 ]; then
# check if there is a more recent year that has already been written
cd ${OUTDIR}
restfile=`ls -t ${ENAM}_restart_Y*.nc | head -1`
suffix="${restfile##*[0-9]}"
number="${restfile%"$suffix"}"
number="${number##*[!-0-9]}"

if [ $YEAR -eq $number ]; then
  echo "Beginning from first year $YEAR"
else
  export YEAR=$number
  echo "Re-initialising from newer restart file '${restfile}' at year $YEAR"
fi
cd ${RUN}
else
echo "Beginning new fresh new run"
fi

# below here is specific for each year:
while [ ${YEAR} -le ${ENDYR} ] ; do
export NEXT=`expr $YEAR + 1`
echo 'YEAR='
echo $YEAR
echo 'NEXT='
echo $NEXT

#source /users/atagliab/xios-2.5s/arch/arch-GCC_BARKLAs.env
source /mnt/data2/users/atagliab/xios-2.5i/arch/arch-GCC_BARKLAifort.env
# clean rundir:
rm -rf ${RUN}/*
#
cp ${BIN}/nemo.exe nemo.exe
cp ${BIN}/namelist_pisces_ref.al.bait namelist_pisces_ref
cp ${BIN}/namelist_pisces_cfg .
cp ${BIN}/namelist_cfg namelist_cfg
cp ${BIN}/namelist_ref .
cp ${BIN}/namelist_top_ref .
cp ${BIN}/namelist_top_cfg_quo_bait namelist_top_cfg
#if (( $YEAR % 10 == 0 )) ;then
if [[ ${YEAR} -eq ${ENDYR} ]];then
cp ${BIN}/file_def_nemo.xml.al.bait file_def_nemo.xml
else
cp ${BIN}/file_def_nemo_EMPTY.xml  file_def_nemo.xml
fi
cp ${BIN}/domain_def_nemo.xml domain_def_nemo.xml
cp ${BIN}/context_nemo.xml context_nemo.xml
cp ${BIN}/field_def_nemo-oce.xml field_def_nemo-oce.xml
cp ${BIN}/field_def_nemo-pisces.xml.al.bait field_def_nemo-pisces.xml
cp ${BIN}/grid_def_nemo.xml grid_def_nemo.xml
cp ${BIN}/iodef.xml iodef.xml
ln -sf ${RUNDIR}/*.nc .
cp ${RUNDIR}/RCP85_CO2_1765_7500.txt .

# can hardcode other ones here if ndeeded:
#
#
# now copy restart
# multi proc restart
#typeset -Z4 i=0
#COUNTER=0
#while [  $COUNTER -le ${NPROC} ]; do
#    i=`printf "%04d" $COUNTER`
#cp ${OUTDIR}/${ENAM}_restart_P${i}_Y${YEAR}.nc restart_trc_${i}.nc
#    let COUNTER=COUNTER+1
#end of processor loop
#done
#
# mono restart
#
cp ${OUTDIR}/${ENAM}_restart_Y${YEAR}.nc restart_trc.nc
# check situatio of rundir
ls -lrt
#######################################
# run model
export OMP_NUM_THREADS=1
#time mpirun -np $XPROC $XBIN/xios_server.exe  : -np $RUNPROC ./nemo.exe
time mpirun -np $RUNPROC ./nemo.exe
# run done
# check if the run was successful (if not, then exit and don't waste my time)
MESSAGE=`tail -40 ocean.output | grep NaN`

if [[ -z $MESSAGE ]]; then
  echo " YEAR $YEAR terminated normally "
else
  echo " E R R O R "
  echo " Abnormal run "
  echo " NaNs in the final statistics of tracers "
  exit
fi
cd ${RUN}
ls -lrt
#source /users/atagliab/xios-2.5s/arch/arch-GCC_BARKLAs.env
source /mnt/data2/users/atagliab/xios-2.5i/arch/arch-GCC_BARKLAifort.env
# copy restart
#./rebuild_nemo PISCES_00001460_restart.trc $RUNPROC
${RBLD2}/rebuild -o PISCES_00001460_restart.trc PISCES_00001460_restart.trc_0???.nc
mv PISCES_00001460_restart.trc.nc ${OUTDIR}/${ENAM}_restart_Y${NEXT}.nc
# save output:
#./rebuild_nemo PISCES_1m_00010101_00011231_ptrc_T $XPROC
#./rebuild_nemo PISCES_1m_00010101_00011231_diad_T $XPROC
#./rebuild_nemo PISCES_1y_00010101_00011231_ptrc_T $XPROC
#./rebuild_nemo PISCES_1y_00010101_00011231_diad_T $XPROC
${RBLD2}/rebuild -o PISCES_1y_00010101_00011231_ptrc_T.nc PISCES_1y_00010101_00011231_ptrc_T_0???.nc
${RBLD2}/rebuild -o PISCES_1m_00010101_00011231_ptrc_T.nc PISCES_1m_00010101_00011231_ptrc_T_0???.nc
${RBLD2}/rebuild -o PISCES_1y_00010101_00011231_diad_T.nc PISCES_1y_00010101_00011231_diad_T_0???.nc
${RBLD2}/rebuild -o PISCES_1m_00010101_00011231_diad_T.nc PISCES_1m_00010101_00011231_diad_T_0???.nc
#${RBLD2}/rebuild -o PISCES_1d_00010101_00011231_bioscaler_T.nc PISCES_1d_00010101_00011231_bioscaler_T_0???.nc
#${RBLD2}/rebuild -o PISCES_1y_00010101_00011231_bioscaler_T.nc PISCES_1y_00010101_00011231_bioscaler_T_0???.nc 
# copy output
mv PISCES_1m_00010101_00011231_ptrc_T.nc ${OUTDIR}/${ENAM}_1m_ptrc_Y${YEAR}.nc
mv PISCES_1m_00010101_00011231_diad_T.nc ${OUTDIR}/${ENAM}_1m_diad_Y${YEAR}.nc
mv PISCES_1y_00010101_00011231_ptrc_T.nc ${OUTDIR}/${ENAM}_1y_ptrc_Y${YEAR}.nc
mv PISCES_1y_00010101_00011231_diad_T.nc ${OUTDIR}/${ENAM}_1y_diad_Y${YEAR}.nc
mv PISCES_1d_00010101_00011231_bioscaler_0001.nc ${OUTDIR}/${ENAM}_1d_bioscaler_Y${YEAR}.nc
mv PISCES_1y_00010101_00011231_bioscalar_0000.nc ${OUTDIR}/${ENAM}_1y_bioscaler_Y${YEAR}.nc
tar -cf ${ENAM}_namelists_Y${YEAR}.tar output.namelist*
mv ${ENAM}_namelists_Y${YEAR}.tar ${OUTDIR}/.
COUNTER=0
#while [  $COUNTER -lt $NPROC ]; do
#    i=`printf "%04d" $COUNTER`
#mv  PISCES_1m_*_ptrc_T_${i}.nc ${OUTDIR}/${ENAM}_1m_ptrc_P${i}_Y${YEAR}.nc
#mv  PISCES_1m_*_diad_T_${i}.nc ${OUTDIR}/${ENAM}_1m_diad_P${i}_Y${YEAR}.nc
#mv  PISCES_1y_*_ptrc_T_${i}.nc ${OUTDIR}/${ENAM}_1y_ptrc_P${i}_Y${YEAR}.nc
#mv  PISCES_1y_*_diad_T_${i}.nc ${OUTDIR}/${ENAM}_1y_diad_P${i}_Y${YEAR}.nc
#mv  PISCES_5d_*_ptrc_T_${i}.nc ${OUTDIR}/${ENAM}_5d_ptrc_P${i}_Y${YEAR}.nc
#mv  PISCES_5d_*_diad_T_${i}.nc ${OUTDIR}/${ENAM}_5d_diad_P${i}_Y${YEAR}.nc
#mv  PISCES_1d_*_ptrc_T_${i}.nc ${OUTDIR}/${ENAM}_1d_ptrc_P${i}_Y${YEAR}.nc
#mv  PISCES_1d_*_diad_T_${i}.nc ${OUTDIR}/${ENAM}_1d_diad_P${i}_Y${YEAR}.nc
#mv  PISCES_*_restart.trc_${i}.nc ${OUTDIR}/${ENAM}_restart_P${i}_Y${NEXT}.nc
#    let COUNTER=COUNTER+1
#done
#
cd ${RUN}
cp ocean.output ${OUTDIR}/${ENAM}_ocean_output_Y${YEAR}
cp co2atm.txt ${OUTDIR}/${ENAM}_co2atm_Y${YEAR}.txt
#end of that particular year
#increment year
export YEAR=$NEXT
#  YEAR=$(( $YEAR + 1 ))
cd ${RUN}
# year loop
echo Job finished date = `date`
done
# create next year script
#  sed s/"export YEAR=3001"/"export YEAR=$ENDYR"/g ${ENAM}.$STYR > ${ENAM}.${ENDYR}a
#  sed s/"export ENDYR=3051"/"export ENDYR=$NEWEND"/g ${ENAM}.${ENDYR}a > ${ENAM}.${ENDYR}
