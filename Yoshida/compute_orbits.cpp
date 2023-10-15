// Assessment 3: n-body gravitational solver

// To avoid warnings tell the compiler to use a recent standard of C++:
// g++ -std=c++17 vector3d.cpp body.cpp compute_orbits.cpp -o compute_orbits
// Run this line to test initially (zero time-steps)
// ./compute_orbits sun_earth.csv test_case.csv 1e-4 0 1

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

#include "vector3d.hpp"
#include "body.hpp"

using std::cout, std::endl;

// *** Function declarations ***

// -----------------------------------------------------------
// ADD YOUR FUNCTION DECLARATIONS HERE

void compute_energy_L(int N, std::vector<body> &system);
void update_acc(int N, int p, std::vector<body> &system);
void yoshida_int(int N, std::vector<body> &system, double dt);

// -----------------------------------------------------------
// Read input data from file
void read_init(std::string input_file, std::vector<body> &system);
// Read the components of a 3d vector from a line
void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z);
// Save the data to file
void save_data(std::ofstream& savefile, const std::vector<body> &system, double t);

int main(int argc, char* argv[])
{
  // Checking if number of arguments is equal to 4:
  if (argc != 6) {
    cout << "ERROR: need 4 arguments - compute_orbits <input_file> <output_file> <dt> <T> <Tsave>" << endl;
    return EXIT_FAILURE;
  }
  // Process command line inputs:
  std::string input_file = argv[1];
  std::string output_file = argv[2];
  double dt = atof(argv[3]); // Time step
  int T = atoi(argv[4]); // Total number of time steps
  int Tsave = atoi(argv[5]); // Number of steps between each save

  std::vector<body> system; // Create an empty vector container for bodies
  read_init(input_file, system); // Read bodies from input file into system
  int N = system.size(); // Number of bodies in system

  cout << "--- Orbital motion simulation ---" << endl;
  cout << " number of bodies N: " << N << endl;
  for(int p=0; p < N; p++)
  {
    cout << "- " << system[p].get_name() << endl; // Display names
  }
  cout << "       time step dt: " << dt << endl;
  cout << "  number of steps T: " << T << endl;
  cout << "   save steps Tsave: " << Tsave << endl;

  std::ofstream savefile (output_file); // Open save file
  if (!savefile.is_open()) {
    cout << "Unable to open file: " << output_file << endl; // Exit if save file has not opened
    return EXIT_FAILURE;
  }
  savefile << std::setprecision(16); // Set the precision of the output to that of a double
  // Write a header for the save file
  savefile << dt << "," << T << "," << Tsave << "," << N << endl;
  for(int p=0; p < (N-1); p++)
  {
    savefile << system[p].get_name() << ",";
  }
  savefile << system[N-1].get_name() << endl;

  // -----------------------------------------------------------
  // ADD YOUR CODE HERE
  
  int count = Tsave;
  double time = 0;
  
    
  for (int t = 0; t < T+1; t++)
  {
      if (count == Tsave)
      {
          compute_energy_L(N, system);
          save_data(savefile, system, time);
//           std::cout << t << std::endl;
          count = 1;
      }
      else {count++;}
      time += dt;
      
      yoshida_int(N, system, dt);
  }

  // -----------------------------------------------------------
  
  savefile.close();
  return EXIT_SUCCESS; 
}

// *** Function implementations ***

// -----------------------------------------------------------
// ADD YOUR FUNCTION IMPLEMENTATIONS HERE

void compute_energy_L(int N, std::vector<body> &system)
{
    double GPE;
    double KE;
    vec L;
    
    // caclulation for the Gravitatinal Potential Energy
    for (int p = 0; p < N; p++)
    {
        KE = 0.5 * system[p].get_mass() * system[p].get_vel().dot(system[p].get_vel());
        L = system[p].get_mass() * system[p].get_pos().cross(system[p].get_vel());
        double GPE = 0;
        
        for (int j = 0; j < N; j++)
        {
            if (j != p)
            {
                GPE -= (system[j].get_mass() * system[p].get_mass()) * (1/( system[j].get_pos() - system[p].get_pos()).length() );
            }
        } 
        
        system[p].set_gpe(GPE);
        system[p].set_ke(KE);
        system[p].set_L(L);
    }
}



void update_acc(int N, int p, std::vector<body> &system)
{

    vec acc(0,0,0); 

    for (int j = 0; j < N; j++)
    {
        if (j != p)
        {
            vec distance = system[j].get_pos() - system[p].get_pos();

            acc += ( system[j].get_mass() / pow(distance.length(), 3) ) * distance;
        }
    }
    system[p].set_acc(acc);
}



void yoshida_int(int N, std::vector<body> &system, double dt)
{
    std::vector<body> new_system = system;
    double w_0 = cbrt(2) / (2 - cbrt(2));
    double w_1 = 1 / (2 - cbrt(2));
    
    
    double c_1 = w_1 / 2;
    double c_2 = (w_0 + w_1) / 2;
    double c_3 = c_2;
    double c_4 = c_1;
    
    double d_1 = w_1;
    double d_2 = w_0;
    double d_3 = d_1;
    
    for (int p = 0; p < N; p++)
    {   
//         1/4
        vec r_i = system[p].get_pos();
        vec v_i = system[p].get_vel();
        
        vec r_1 = r_i + (c_1 * v_i * dt);
        new_system[p].set_pos(r_1);
        
        update_acc(N, p, new_system);
        vec a_1 = new_system[p].get_acc();
        vec v_1 = v_i + (d_1 * a_1 * dt);
        
//         2/4
        vec r_2 = r_1 + (c_2 * v_1 * dt);
        new_system[p].set_pos(r_2);
        
        update_acc(N, p, new_system);
        vec a_2 = new_system[p].get_acc();
        vec v_2 = v_1 + (d_2 * a_2 * dt);
        
//         3/4
        vec r_3 = r_2 + (c_3 * v_2 * dt);
        new_system[p].set_pos(r_3);
        
        update_acc(N, p, new_system);
        vec a_3 = new_system[p].get_acc();
        vec v_3 = v_2 + (d_3 * a_3 * dt);
        
//         4/4
        vec r_4 = r_3 + (c_4 * v_3 * dt);
        vec v_4 = v_3;
        
        system[p].set_pos(r_4);
        system[p].set_vel(v_4);
    }
}



// -----------------------------------------------------------
void read_init(std::string input_file, std::vector<body> &system)
{
  std::string line; // Declare a string to store each line
  std::string name; // String to store body name
  double m, x, y, z, vx, vy, vz; // Doubles to store vector components
  int line_cnt = 0; // Line counter

  // Declare and initialise an input file stream object
  std::ifstream data_file(input_file); 

  while (getline(data_file, line)) // Read the file line by line
  {
    line_cnt++;
    std::stringstream data_line(line); // Create a string stream from the line
    switch (line_cnt)
    {
      case 1:
        name = line;
        break;
      case 2:
        m = std::stod(line); // Convert string line into double
        break;            
      case 3:
        read_vector3d(data_line, x, y, z); // Read the 3 components of the vector on the line
        break;
      case 4:
        read_vector3d(data_line, vx, vy, vz); // Read the 3 components of the vector on the line
        break;                   
      }
      if (line_cnt==4) // Data for one body has been extracted
      {
        line_cnt = 0; // Reset line counter
        body b(name,m,vec(x,y,z),vec(vx,vy,vz)); // Package data into body
        system.push_back(b); // Add body to system
      }
    }
    // Close the file
    data_file.close();
}

void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z)
{
  std::string value; // Declare a string to store values in a line
  int val_cnt = 0; // Value counter along the line
  
  while (getline(data_line, value, ','))
  {
    val_cnt++;
    switch (val_cnt)
    {
      case 1:
        x = std::stod(value); 
        break;
      case 2:
        y = std::stod(value); 
        break;
      case 3:
        z = std::stod(value);
        break;
    }
  } 
}

void save_data(std::ofstream& savefile, const std::vector<body> &system, double t)
{
    // Function for saving the simulation data to file.

    vec L; // Total angular momentum
    double E = 0.0; // Total energy
    for(int p = 0; p < system.size(); p++)
    { 
      E += system[p].get_ke() + 0.5*system[p].get_gpe();
      L += system[p].get_L();
    }
    double Lmag = L.length(); // Magnitude of total angular momentum

    // Write a header for this time-step with time, total energy and total mag of L:
    savefile << t << "," << E << "," << Lmag << endl;

    // Loop over the bodies:
    for(int p = 0; p < system.size(); p++)
    { 
      // Output position and velocity for each body:
      savefile << system[p].get_pos() << "," << system[p].get_vel() << endl;
    }
}