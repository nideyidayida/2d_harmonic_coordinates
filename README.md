# Final Project: Harmonic Coordinates for Character Animation

## Clone the repository:
git clone https://github.com/NYUGP17/final_project_jockyxu1122  
## Build, compile, and run the project:
cd final_project_jockyxu1122  
mkdir build  
cd build  
cmake -DCMAKE_BUILD_TYPE=Release ../  
make  
./final_project_bin *some_2D_mesh_file*

## UI:
Enter SELECT mode to pin cage points.  
After seleting all cage points, press 'H'.  
Enter TRANSLATE mode then drag any cage vertex.
### Key:
'S': enter SELECT mode.  
'H': compute H function.  
'M': enter TRANSLATE mode.  
'N': enter NONE mode.
