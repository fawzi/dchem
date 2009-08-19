/// the recommended fundamental constants of physics and chemistry
/// based on the CODATA 2006
///
///      Fundamental physical constants (SI units)
///      Literature: - B. J. Mohr and B. N. Taylor,
///                    "CODATA recommended values of the fundamental physical
///                     constants: 1998", Rev. Mod. Phys. 72(2), 351 (2000).
///                  - B. J. Mohr and B. N. Taylor,
///                    "CODATA recommended values of the fundamental physical
///                     constants: 2002", Rev. Mod. Phys. 77(2), 1 (2005).
///                  - B. J. Mohr and B. N. Taylor,
///                    "CODATA recommended values of the fundamental physical
///                     constants: 2006", http://physics.nist.gov/constants
///
/// based on the cp2k version by JGH, MK & AK Copyright (C) 2000 - 2009  CP2K developers group
/// fawzi
/// License: GPLv3
module dchem.Physcon;
public import tango.math.Math:pi;

static this(){
/// Planck constant [J*s]
const real h_planck = 6.62606896E-34;
const real h_bar = h_planck/(2.0*pi)

/// Elementary charge [C]
const real e_charge = 1.602176487E-19;

/// Electron mass [kg]
const real e_mass = 9.10938215E-31;

/// Proton mass [kg]
const real p_mass = 1.672621637E-27;

/// Electron g factor [ ]
const real e_gfactor = -2.0023193043622;

/// Fine-structure constant
//MK a_fine = 0.5*mu_perm*c_light*e_charge**2/h_planck
const real a_fine = 7.2973525376E-3;

/// Rydberg constant [1/m]
//MK rydberg = 0.5*e_mass*c_light*a_fine**2/h_planck
const real rydberg = 10973731.568527;

/// Avogradro constant [1/mol]
const real n_avogadro = 6.02214179E+23;

/// Boltzmann constant [J/K]
const real boltzmann = 1.3806504E-23;

/// Atomic mass unit [kg]; conversion factor [u] -> [kg]
const real a_mass = 1.660538782E-27;

/// Bohr radius [m]
//MK a_bohr = a_fine/(4.0*pi*rydberg)
const real a_bohr = 0.52917720859E-10;

// Conversion factors

/// [u] -> [a.u.]
const real massunit = a_mass/e_mass;

/// [Bohr] -> [Angstrom]
const real angstrom = 1.0E+10*a_bohr;

/// [Angstrom] -> [Bohr]
const real bohr = 1.0/angstrom;

/// [a.u.] -> [s]
const real seconds = 1.0/(4.0*pi*rydberg*c_light);

/// [a.u.] -> [fs]
const real femtoseconds = 1.0E+15*seconds;

/// [a.u.] -> [ps]
const real picoseconds = 1.0E+12*seconds;

/// [a.u.] -> [J]
const real joule = 2.0*rydberg*h_planck*c_light;

/// [a.u.] -> [K]
const real kelvin = joule/boltzmann;

/// [a.u.] -> [kJ/mol]
const real kjmol = 0.001*joule*n_avogadro;

/// [a.u.] -> [kcal/mol]
const real kcalmol = kjmol/4.184;

/// [a.u.] -> [Pa]
const real pascal = joule/a_bohr**3;

/// [a.u.] -> [bar]
const real bar = pascal/1.0E+5;

/// [a.u.] -> [atm]
const real atm = pascal/1.013250E+5;

/// [a.u.] -> [eV]
const real evolt = joule/e_charge;

/// [a.u.] -> [Hz]
const real hertz = joule/h_planck;

// [a.u./Bohr**2] -> [1/cm] (wave numbers)
const real vibfac;
static this(){
    vibfac = 5.0*sqrt(kjmol)/(pi*a_bohr*c_light);
}

/// [a.u.] -> [1/cm] (wave numbers)
const real wavenumbers = 0.02*rydberg;

/// [a.u.] -> [esu] (electrostatic units)
const real[3] esu;
/// [a.u.] -> [debye] (electrostatic units)
const real debye;
static this(){
    esu[0] = 1.0E+21*a_bohr*c_light*e_charge;
    for (int i=1;i<esu.length;++i){
        esu[i] = esu[i-1]/bohr;
    }
    debye = esu(1);
}

