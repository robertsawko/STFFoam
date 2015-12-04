#!/bin/bash
#
#PBS -N RSbubble
#PBS -k oe
#PBS -l select=2:ncpus=16
#PBS -q medium
#PBS -m abe -M r.sawko@cranfield.ac.uk

#Load OF environment
#OpenFOAM can be loaded in different ways: single/ double 
#precision or with different MPI implementations. This 
#script covers INTELMPI and SGIMPI. INTELMPI is somewhat
#simpler to run but here it commented out.
#. /usr/local/non-commercial/OpenFOAM/OpenFOAM-SGISLES11-2.1.X_X8664_21JUNE12_REPO/setup.sh DP INTELMPI
. /usr/local/non-commercial/OpenFOAM/OpenFOAM-SGISLES11-2.1.X_X8664_21JUNE12_REPO/setup.sh DP SGIMPI

#Creating some variables for SGIMPI run.
#Not all the variables here are necessary.
cpus=$( cat $PBS_NODEFILE | wc -l)
FOAM_NODES=$(cat $PBS_NODEFILE )
FIRSTNODE=`echo $FOAM_NODES|awk '{print $1}'`
FOAM_NODES=$( echo $FOAM_NODES | sed -e 's/ /, /g')
MPI_UNIVERSE=`echo "$FOAM_NODES $NCPUS"`
export MPI_UNIVERSE=$MPI_UNIVERSE
TOTCPUS=`echo $NCPUS*$cpus|bc`

#Echo some information on acquired compute nodes
echo "Number of compute nodes: " $cpus 
echo "Names of acquired compute nodes: " $FOAM_NODES
echo "Total number of cores used by this job: " $TOTCPUS

cd $PBS_O_WORKDIR
##########Start your OpenFOAM scripting here##########
#Put all the applications that need to be run prior 
#to the parallel solver run. This way you copy only
#case setting files between computers.
#e.g.
#rm -r processor*
cp -r 0.org 0
#Create mesh
blockMesh > log.blockMesh 2> err.blockMesh
#Decompose the domain
decomposePar > log.decomposePar 2> err.decomposePar

#####################Parallel run#####################
#IntelMPI or SGIMPI

#For IntelMPI hostfile switch works but you need to 
#calculate the total number of nodes. It might be better
#to set it up manually!

#mpirun -hostfile $PBS_NODEFILE -np $TOTCPUS simpleFoam -parallel 

#For SGIMPI you need to set MPI_UNIVERSE. By default mpirun on the
#localhost.
mpirun $MPI_UNIVERSE TFRInterFoam -parallel > log.TFRInterFoam 2>err.TFRInterFoam
