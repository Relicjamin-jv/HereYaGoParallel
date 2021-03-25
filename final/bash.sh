#Runs the bash for the serial and MPI and also outs that output to a file

echo -e "Collin Campbell ~ COMP 233A ~ serial and MPI Bash Script\n" > timingInfo.txt
echo -e "Running the serial"
./jacobiSerial .01 100000 serialPPM.ppm 
echo -e "Running 4 processes"
mpirun -np 4 ./MPIJacobi .01 100000 mpiPPM4.ppm  #runs the MPI code 
echo -e "Running 8 processes"
mpirun -np 8 ./MPIJacobi .01 100000 mpiPPM8.ppm  #runs the MPI code 
echo -e "Running 16 processes"
mpirun -np 16 ./MPIJacobi .01 100000 mpiPPM16.ppm  #runs the MPI code 
echo -e "\n\nBash Script Terminated\n" >> timingInfo.txt
