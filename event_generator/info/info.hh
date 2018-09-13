#ifndef _INFO_HH_
#define _INFO_HH_

const double amu = 931.494;// MeV
const int		 MAX_STORY = 50000; 	 // Total number of generated events
//const double ENERGY 	 = 397.6827; // Beam energy (MeV/u)
const double ENERGY 	 = 200; // Beam energy (MeV/u)
// const int		 A 				 = 12;			 // Mass number of the incident nucleus A
const int		 A 				 = 11;			 // Mass number of the incident nucleus A
const double Exe  		 = 5.00000;  // Residual excitation energy (MeV) (change it for deeply bound states)
//const double MOM_SIGMA = 105.00;   // Internal momentum spread (Gauss)
const double MOM_SIGMA = 50.00;   // Internal momentum spread (Gauss)
const bool   ISOTROPIC = true;   	 // proton-proton CM scattering

// const double MA  = 11174.950;			 // Nuclear mass of incident A nucleus(MeV/c²)
const double MA  = 11.0421*amu;			 // Nuclear mass of incident A nucleus(MeV/c²)
// const double MB  = 10252.628 + Exe;// Nuclear mass of the residual fragment B (MeV/c²)
const double MBgs  = 10.0517*amu;// Nuclear mass of the residual fragment B (MeV/c²)
const double MB  = 10.0517*amu + Exe;// Nuclear mass of the residual fragment B (MeV/c²)
const double Ma  = 938.279;  			 // Mass of the knocked-out nucleon (MeV/c²)
const double Mi  = 938.279;  			 // Mass of the scattered (target) nucleon (MeV/c²)

//Constants
const double PI  = 3.14159265358979323846;

#endif
