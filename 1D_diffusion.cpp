#include "base.h"

void writeVTK(  MultiBlockLattice2D<T,ADESCRIPTOR>& adLattice,
                SimulationParams<T> & simParams,
                plint iter,
				string f_name)
{
    //T dx = parameters.getDeltaX();
    T dx = 1.0/simParams.getResolution();
	plint resolution = simParams.getResolution();
	unique_ptr<MultiScalarField2D<T> > density = computeDensity(adLattice);
	unique_ptr<MultiTensorField2D<T,2> > gradient = computeGradient(*density);
	//unique_ptr<MultiScalarField2D<T> > flux = computeMassFlux(adLattice,boundary_domain,simParams);
    VtkImageOutput2D<T> vtkOut(createFileName(f_name.c_str(), iter, 6), dx);
    vtkOut.writeData<float>(*density, "density", 1.);
    vtkOut.writeData<2,float>(*gradient, "potential flux",resolution);
    vtkOut.writeData<float>(*computeNorm(*gradient), "norm flux",resolution);
    vtkOut.writeData<2,float>(*computeGradient(*computeNorm(*gradient)), "Laplacian",resolution);
    pcout << "Writing vtk_instance at iT = " << iter << endl;
    pcout << "Average Flux = " << computeAverage(*computeGradientNorm(*density), adLattice.getBoundingBox()) << endl;
}


void domain_setup(  MultiBlockLattice2D<T,ADESCRIPTOR>& adLattice,
                    SimulationParams<T> & simParams)
{
    // Initial condition and constant concentration boundary condition at south
    // boundary 

    plint nx = simParams.getNx();
    T c_init = 0.0;
    Box2D south_boundary(0,nx-1,0,0);

    initializeAtEquilibrium(    adLattice,
                                adLattice.getBoundingBox(),
                                c_init,
                                Array<T,2>((T)0.0,(T)0.0));

    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, ADESCRIPTOR>* 
        bc = createLocalAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>();

    bc->addTemperatureBoundary1N(south_boundary,adLattice);

    setBoundaryDensity(adLattice,south_boundary,(T)1.0);

    adLattice.initialize();
}


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    // Instantiate variables
    SimulationParams<T> simParams = assign_params("params.xml");

    T Omega = 1.0/simParams.getTau_ad();
    pcout << "Omega = " << Omega << endl;
    pcout << "tau_ad (tau advection/diffusion lattice) = " << simParams.getTau_ad() << endl;
    plint nx = simParams.getNx();
    plint ny = simParams.getNy();
    plint resolution = simParams.getResolution();

    // Instantiate Charge Carrier Lattice
    MultiBlockLattice2D<T,ADESCRIPTOR> adLattice(
                                    nx,
                                    ny,
                                    new ADYNAMICS<T, ADESCRIPTOR>(Omega));

    domain_setup(  adLattice, simParams);

    plint convergenceIter = simParams.getConvergenceIter();
    // =========================================================================
    // Time Stepping
    // =========================================================================
    for (plint iT = 0; iT<=simParams.getMaxIter(); iT++)
    {
        if (iT%convergenceIter==0) {
            pcout << "Writing vtk_instance at iT = " << iT << endl;
            writeVTK(   adLattice, simParams, iT, "Ad_Diffusion Lattice");
        }
        adLattice.collideAndStream();
    }
    plb_ofstream succesiveProfiles("concentration_final.dat");
    succesiveProfiles << std::setprecision(7)
        << *computeDensity(adLattice,Box2D(0,0,0,resolution)) << endl <<endl;
}