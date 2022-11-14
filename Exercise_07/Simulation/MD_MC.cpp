#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "MD_MC.h"

using namespace std;

int main()
{ 

//SOLID SIMULATION
  string input = "input_solid.in";
  string config = "config.fcc";
  string output_solid[2] = {"../Solid_sim/output_epot.dat", "../Solid_sim/output_press.dat"};
  string output_radial = "../Solid_sim/output_radial.dat";
  string output_radial_Verlet = "../Solid_sim/output_radial_Verlet.dat";
  Input(input, config); //Inizialization
  string out_eq = "../Solid_sim/Equi_solid.dat";
  ofstream out_ac_err;
  out_ac_err.open("../Solid_sim/Autocorr_err_solid.dat");
  int equi_time_solid = 500;
  Energy_values_with_print(out_eq, equi_time_solid); 
  //Inst energy values are printed in order to understand the needed equilibration time
  //int nconf = 1;
  Input(input, config);
  //Equilibration
  for(int j = 0; j < 200; j++){
    Move();
  }
  //Inst energy values are printed in order to understand the autocorrelation function and the error trend 
  //with data block sizing

  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure();
      out_ac_err << nstep*(iblk-1) + istep << "\t" << walker[iv]/(double)npart << endl;
      Measure_radial();
      Accumulate(); //Update block averages
    }
    Averages_Metropolis(iblk, output_solid);   //Print results for current block

    Average_radial(iblk, output_radial);
  }
  input = "input_solid_Verlet.in";
  Input(input, config); //Inizialization for Verlet algorithm

  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    if(iblk == 1){
      for(int j = 0; j < 2000; j++){
        Move();
      }
    }
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure_radial();
      Accumulate(); //Update block averages
    }

    Average_radial(iblk, output_radial_Verlet);
  }
  out_ac_err.close();
//END OF SOLID SIMULATION

//LIQUID SIMULATION
  input = "input_liquid.in";
  config = "config.fcc";
  string output_liquid[2] = {"../Liquid_sim/output_epot.dat", "../Liquid_sim/output_press.dat"};
  output_radial = "../Liquid_sim/output_radial.dat";
  output_radial_Verlet = "../Liquid_sim/output_radial_Verlet.dat";
  Input(input, config); //Inizialization
  out_eq = "../Liquid_sim/Equi_liquid.dat";
  out_ac_err.open("../Liquid_sim/Autocorr_err_liquid.dat");
  int equi_time_liquid = 500;
  Energy_values_with_print(out_eq, equi_time_liquid);
  //Inst energy values are printed in order to understand autocorr function and trend of error qith variou data-blocks size
  //int nconf = 1;
  Input(input, config);
  for(int j = 0; j < 200; j++){
    Move();
  } 
  //Equilibration  
  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure();
      out_ac_err << nstep*(iblk-1) + istep << "\t" << walker[iv]/(double)npart << endl;
      Measure_radial();
      Accumulate(); //Update block averages
    }
    Averages_Metropolis(iblk, output_liquid);   //Print results for current block
    Average_radial(iblk, output_radial);
  }

  input = "input_liquid_Verlet.in";
  Input(input, config); //Inizialization for Verlet algorithm

  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    if(iblk == 1){
      for(int j = 0; j < 2000; j++){
        Move();
      }
    }
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure_radial();
      Accumulate(); //Update block averages
    }

    Average_radial(iblk, output_radial_Verlet);
  }
  out_ac_err.close();
//END OF LIQUID SIMULATION

//GAS SIMULATION
  input = "input_gas.in";
  config = "config.fcc";
  string output_gas[2] = {"../Gas_sim/output_epot.dat", "../Gas_sim/output_press.dat"};
  output_radial = "../Gas_sim/output_radial.dat";
  output_radial_Verlet = "../Gas_sim/output_radial_Verlet.dat";
  Input(input, config); //Inizialization
  out_eq = "../Gas_sim/Equi_gas.dat";
  out_ac_err.open("../Gas_sim/Autocorr_err_gas.dat");
  int equi_time_gas = 500;
  Energy_values_with_print(out_eq, equi_time_gas);

  Input(input, config);
  for(int j = 0; j < 200; j++){
    Move();
  } 
  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure();
      out_ac_err << nstep*(iblk-1) + istep << "\t" << walker[iv]/(double)npart << endl;
      Measure_radial();
      Accumulate(); //Update block averages
    }
    Averages_Metropolis(iblk, output_gas);   //Print results for current block
    Average_radial(iblk, output_radial);
  }

  input = "input_gas_Verlet.in";
  Input(input, config); //Inizialization for Verlet algorithm

  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    if(iblk == 1){
      for(int j = 0; j < 5000; j++){
        Move();
      }
    }
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure_radial();
      Accumulate(); //Update block averages
    }

    Average_radial(iblk, output_radial_Verlet);
  }
  out_ac_err.close();
//END OF GAS SIMULATION*/

  ConfFinal(); //Write final configuration

  return 0;
}


void Energy_values_with_print(string filename, int time){
  ofstream out;
  out.open(filename);
  for (int i = 0; i < time; i++){
    Move();
    Measure_energy();
    out << i << "\t" << walker[iv]/(double)npart<< endl;
  }
  out.close();
}

void Equilibration(int time){
  for (int i = 0; i < time; i++){
  Move();
  }
}

void Input(string input, string config)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open("../../Random/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open(input);

  ReadInput >> iNVET;
  ReadInput >> restart;

  if(restart) Seed.open("seed.out");
  else Seed.open("../../Random/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  nbins = m_props-n_props;
  bin_size = box/(2*nbins);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  iw = 4; //Pressure
  n_props = 5; //Number of observables

  for(int i = 0; i < m_props; i++){
    walker[i] = 0;
  }

//Read initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart)
  {
    ReadConf.open("config.out");
    ReadVelocity.open("velocity.out");
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else 
  {
    ReadConf.open(config);
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
//Evaluate properties of the initial configuration
  Measure();

//Print initial values for measured properties
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  if(iNVET == 0) cout << "Initial temperature      = " << walker[it] << endl;
  if(iNVET == 0) cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  if(iNVET == 0) cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Initial pressure         = " << walker[iw]/(double)npart << endl;

  return;
}

void Move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= or <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement
{
  double v = 0.0, w = 0.0, kin=0.0;
  double vij, wij;
  double dx, dy, dz, dr;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
        v += vij;
        w += wij;
      }
    }          
  }


  for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  double pot_tail = 8*M_PI*rho/(3.*pow(rcut,3))*(1./(3.*pow(rcut,6))-1);
  double press_tail = 32*M_PI*rho/(3.*pow(rcut,3))*(1./(3*pow(rcut,6))-1./2);
  walker[iv] = 4.0 * v + pot_tail; // Potential energy
  if(iNVET == 0) walker[ik] = kin; // Kinetic energy
  if(iNVET == 0) walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  if(iNVET == 0) walker[ie] = 4.0 * v + kin;  // Total energy
  walker[iw] = 16.0 * w/vol + rho * temp + press_tail; //Pressure

  return;
}

void Measure_energy() //Energy meas
{
  double v = 0.0;
  double vij;
  double dx, dy, dz, dr;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        v += vij;
      }
    }          
  }


  double pot_tail = 8*M_PI*rho/(3.*pow(rcut,3))*(1./(3.*pow(rcut,6))-1);
  walker[iv] = 4.0 * v + pot_tail; // Potential energy

  return;
}

void Measure_radial() //measure radial distrib function
{
  double dx, dy, dz, dr;

  for(int j = n_props; j < m_props; j++){
    walker[j] = 0;
  }

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < box/2.){
        int k = int(dr/bin_size) + n_props;
        walker[k] += 2;
      }
    }          
  }

  return;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<m_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<m_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<m_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk, string file_out[5]) //Print results for current block
{
    
   ofstream Epot, Ekin, Etot, Temp, Press;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open(file_out[iv],ios::app);
    Ekin.open(file_out[ik],ios::app);
    Temp.open(file_out[it],ios::app);
    Etot.open(file_out[ie],ios::app);
    Press.open(file_out[iw],ios::app);
    
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp = Error(glob_av[it],glob_av2[it],iblk);

    stima_press = blk_av[iw]/blk_norm; //Pressure
    glob_av[iw] += stima_press;
    glob_av2[iw] += stima_press*stima_press;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);



//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Kinetic energy
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Total energy
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
//Pressure
    Press << setw(wd) << iblk <<  setw(wd) << stima_press << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Press.close();

}

void Averages_Metropolis(int iblk, string file_out[2]) //Print results for current block
{
    
   ofstream Epot, Press;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
  
    Epot.open(file_out[0],ios::app);
    Press.open(file_out[1],ios::app);

    
    stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_press = blk_av[iw]/blk_norm; //Pressure
    glob_av[iw] += stima_press;
    glob_av2[iw] += stima_press*stima_press;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);



//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Pressure
    Press << setw(wd) << iblk <<  setw(wd) << stima_press << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Press.close();

}

void Average_radial(int iblk, string file_out) //Print results for current block
{
    
    ofstream Rad;
    const int wd=12;
    
    if (iNVET == 0) cout << "Block number " << iblk << endl;
    if (iNVET == 0) cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Rad.open(file_out,ios::app);
    for(int h = n_props; h < m_props; h++){
      double del_V = 4*M_PI/3.*(pow ((h-n_props+1)*bin_size ,3) - pow ((h-n_props)*bin_size ,3));
      stima_rad = blk_av[h]/(rho*m_part*del_V)/blk_norm; //radial
      glob_av[h] += stima_rad;
      glob_av2[h] += stima_rad*stima_rad;
      err_rad=Error(glob_av[h],glob_av2[h],iblk);
      if(iblk == nblk){
        
        Rad << setw(wd) << iblk <<  setw(wd) << stima_rad << setw(wd) << glob_av[h]/(double)iblk << setw(wd) << err_rad << endl;
      }
    }



    if (iNVET == 0) cout << "----------------------------" << endl << endl;

    Rad.close();

}


void ConfFinal(void)
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteVelocity.open("velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}
