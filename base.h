#include "palabos2D.h"
#include "palabos2D.hh"
#include <cmath>

using namespace plb;
using namespace std;

typedef double T;

#define ADESCRIPTOR descriptors::AdvectionDiffusionD2Q5Descriptor
#define ADYNAMICS AdvectionDiffusionBGKdynamics

// Gas Constant and Faraday's Constant
T const Ru = 8.3144598;
T const F = 96485.33289;

//------------------------------------------------------------------------------
//  These are debugging functions for outputting they are separate from IO so
//  that IO can be defined downstream
//------------------------------------------------------------------------------
// {{{ printPopulations
void printPopulations(	plint iX,
						plint iY,
						BlockLattice2D<T,ADESCRIPTOR>& adLattice,
						string popName)
{
	Array<T,5>& cell = adLattice.get(iX,iY).getRawPopulations();
	pcout << "printing " << popName << " @ ix,iy = " << iX << "," << iY << endl;
	T density = 0.0;
	for (int i = 0; i < 5; i++)
	{
		pcout << "i = " << i << "; population = " << cell[i]+ADESCRIPTOR<T>::t[i] << endl;
		density+=cell[i]+ADESCRIPTOR<T>::t[i];
	}
	pcout << "concentration =  " << density << endl;
}
// }}}
// {{{  print_array 
template<int nrows, int ncols>
void print_array(   T const X[nrows][ncols],
                    string arr_name){
    pcout << arr_name << " = " << endl;
    for (int i = 0; i<nrows; i++){
        for (int j = 0; j<ncols; j++){
            pcout << "--> " << arr_name << "[" << i << "][" << j << "] = " << X[i][j] << endl;
        }
    }
}
// }}}


//------------------------------------------------------------------------------
// Essentially All the methods need SimulationParams so it is defined at the top
// of the stream as well as its I/O functions
//------------------------------------------------------------------------------
// {{{ start SimulationParams
template<typename T>
class SimulationParams
{
	private:
		plint lx, ly, resolution, convergenceIter, nx, ny, maxIter;
		T tau_ad;
	public:
		SimulationParams<T>(
					plint lx_,
					plint ly_,
					plint resolution_,
					plint convergenceIter_,
					T tau_ad_
					):
			lx(lx_),
			ly(ly_),
			resolution(resolution_),
			convergenceIter(convergenceIter_),
			tau_ad(tau_ad_)
		{
			//Calculated values
			nx = lx*resolution+1;
			ny = ly*resolution+1;
			maxIter = (resolution)*(resolution);
		}
		//Methods
		plint	getLx() const { return lx; }
		plint	getLy() const { return ly; }
		plint	getResolution() const { return resolution; }
		plint	getConvergenceIter() const { return convergenceIter; }
		plint	getNx() const { return nx; }
		plint	getNy() const { return ny; }
		plint	getMaxIter() const { return maxIter; }
		T		getTau_ad() const { return tau_ad; }
};
// }}} end SimulationParams
// {{{ start assign_params
SimulationParams<T> assign_params(string f_name)
{
	plint lx, ly, resolution, convergenceIter;
	T tau_ad;
	try{
		XMLreader xmlFile(f_name);
		pcout << "CONFIGURATION" << endl;
		pcout << "=============" << endl;
		xmlFile.print(0);
		xmlFile["inputs"]["lx"].read(lx);
		xmlFile["inputs"]["ly"].read(ly);
		xmlFile["inputs"]["resolution"].read(resolution);
		xmlFile["inputs"]["convergenceIter"].read(convergenceIter);
		xmlFile["inputs"]["tau_ad"].read(tau_ad);
		pcout << "=============" << endl << endl;
		} catch (PlbIOException& exception) { 
		    pcout << exception.what() << endl;
		}
		return SimulationParams<T> (
					lx,
					ly,
					resolution,
					convergenceIter,
					tau_ad);
};
//}}} end assign_params
// {{{ start writeLogFile
template<typename T>
void writeLogFile(  IncomprFlowParam<T> const& parameters,
		SimulationParams<T> const& simParams,
		std::string const& title)
{
	std::string fullName = global::directories().getLogOutDir() + "plbLog.dat";
	plb_ofstream ofile(fullName.c_str());
	ofile << title << "\n\n";
	ofile << "Velocity in lattice units: u=" << parameters.getLatticeU() << "\n";
	ofile << "Reynolds number:           Re=" << parameters.getRe() << "\n";
	ofile << "Lattice resolution:        N=" << parameters.getResolution() << "\n";
	ofile << "Relaxation frequency:      omega=" << parameters.getOmega() << "\n";
	ofile << "LatticeViscosity:          nu=" << parameters.getLatticeNu() << "\n";
	ofile << "Extent of the system:      lx=" << parameters.getLx() << "\n";
	ofile << "nx:                        nx=" << parameters.getNx() << "\n";
	ofile << "ny:                        ny=" << parameters.getNy() << "\n";
	ofile << "Extent of the system:      ly=" << parameters.getLy() << "\n";
	ofile << "Extent of the system:      lz=" << parameters.getLz() << "\n";
	ofile << "Grid spacing deltaX:       dx=" << parameters.getDeltaX() << "\n";
	ofile << "Time step deltaT:          dt=" << parameters.getDeltaT() << "\n";
	ofile << "==============================" << "\n";
	ofile << "lx			lx = " << simParams.getLx() << "\n";
	ofile << "ly			ly = " << simParams.getLy() << "\n";
	ofile << "resolution			resolution = " << simParams.getResolution() << "\n";
	ofile << "convergenceIter			convergenceIter = " << simParams.getConvergenceIter() << "\n";
	ofile << "nx			nx = " << simParams.getNx() << "\n";
	ofile << "ny			ny = " << simParams.getNy() << "\n";
	ofile << "maxIter			maxIter = " << simParams.getMaxIter() << "\n";
	ofile << "tau_ad			tau_ad = " << simParams.getTau_ad() << "\n";
}
// }}} end writeLogFile


















