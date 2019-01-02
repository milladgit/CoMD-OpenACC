/// \file
/// Leapfrog time integrator

#include "timestep.h"

#include <omp.h>

#include "CoMDTypes.h"
#include "linkCells.h"
#include "parallel.h"
#include "performanceTimers.h"

#include "accFunctions.h"
#include <openacc.h>

static void advanceVelocity(SimFlat* s, int nBoxes, real_t dt);
static void advancePosition(SimFlat* s, int nBoxes, real_t dt);


/// Advance the simulation time to t+dt using a leap frog method
/// (equivalent to velocity verlet).
///
/// Forces must be computed before calling the integrator the first time.
///
///  - Advance velocities half time step using forces
///  - Advance positions full time step using velocities
///  - Update link cells and exchange remote particles
///  - Compute forces
///  - Update velocities half time step using forces
///
/// This leaves positions, velocities, and forces at t+dt, with the
/// forces ready to perform the half step velocity update at the top of
/// the next call.
///
/// After nSteps the kinetic energy is computed for diagnostic output.
double timestep(SimFlat* s, int nSteps, real_t dt)
{
	for (int ii=0; ii<nSteps; ++ii)
	{
		startTimer(velocityTimer);
		advanceVelocity(s, s->boxes->nLocalBoxes, 0.5*dt); 
		stopTimer(velocityTimer);

		startTimer(positionTimer);
		advancePosition(s, s->boxes->nLocalBoxes, dt);
		stopTimer(positionTimer);

		startTimer(redistributeTimer);
		redistributeAtoms(s);
		stopTimer(redistributeTimer);

		startTimer(computeForceTimer);
		computeForce(s);
		stopTimer(computeForceTimer);

		startTimer(velocityTimer);
		advanceVelocity(s, s->boxes->nLocalBoxes, 0.5*dt); 
		stopTimer(velocityTimer);
	}

	kineticEnergy(s, 1);

	return s->ePotential;
}

void computeForce(SimFlat* s)
{
	s->pot->force(s);
}


void advanceVelocity(SimFlat* s, int nBoxes, real_t dt)
{
	#ifdef __pointerchain
	#pragma pointerchain declare(s->species{SpeciesData*}, s->boxes->nTotalBoxes{int}, s->boxes->nbrBoxes{int**})
	#pragma pointerchain declare(s->atoms->r{real3*}, s->atoms->p{real3*}, s->atoms->f{real3*}, s->atoms->U{real_t*})
	#pragma pointerchain declare(s->boxes->nAtoms{int*}, s->atoms->iSpecies{int*}, s->atoms->gid{int*})
	#endif

	int fSize = s->boxes->nTotalBoxes*MAXATOMS;

	#ifdef __pointerchain
	#pragma pointerchain region begin
	#endif
	#pragma acc parallel loop gang independent present(s->boxes->nAtoms[0:s->boxes->nTotalBoxes], \
		s->atoms->p[0:fSize][0:3], s->atoms->f[0:fSize][0:3])
	for (int iBox=0; iBox<nBoxes; iBox++)
	{
		#pragma acc loop vector independent
		for (int ii=0; ii<s->boxes->nAtoms[iBox]; ii++)
		{
			int iOff = MAXATOMS * iBox + ii;

			s->atoms->p[iOff][0] += dt*s->atoms->f[iOff][0];
			s->atoms->p[iOff][1] += dt*s->atoms->f[iOff][1];
			s->atoms->p[iOff][2] += dt*s->atoms->f[iOff][2];
		}
	}
	#ifdef __pointerchain
	#pragma pointerchain region end
	#endif
}

void advancePosition(SimFlat* s, int nBoxes, real_t dt)
{
	#ifdef __pointerchain
	#pragma pointerchain declare(s->species{SpeciesData*}, s->boxes->nTotalBoxes{int}, s->boxes->nbrBoxes{int**})
	#pragma pointerchain declare(s->atoms->r{real3*}, s->atoms->p{real3*}, s->atoms->f{real3*}, s->atoms->U{real_t*})
	#pragma pointerchain declare(s->boxes->nAtoms{int*}, s->atoms->iSpecies{int*}, s->atoms->gid{int*})
	#endif

	int fSize = s->boxes->nTotalBoxes*MAXATOMS;

	#ifdef __pointerchain
	#pragma pointerchain region begin
	#endif
	#pragma acc parallel loop gang independent present(s->boxes->nAtoms[0:s->boxes->nTotalBoxes], \
		s->atoms->p[0:fSize][0:3], s->atoms->r[0:fSize][0:3], s->atoms->iSpecies[0:fSize], \
		s->species[0:1])
	for (int iBox=0; iBox<nBoxes; iBox++)
	{
		#pragma acc loop vector independent
		for (int ii=0; ii<s->boxes->nAtoms[iBox]; ii++)
		{
			int iOff = MAXATOMS * iBox + ii;
			int iSpecies = s->atoms->iSpecies[iOff];
			real_t invMass = 1.0/s->species[iSpecies].mass;
			s->atoms->r[iOff][0] += dt*s->atoms->p[iOff][0]*invMass;
			s->atoms->r[iOff][1] += dt*s->atoms->p[iOff][1]*invMass;
			s->atoms->r[iOff][2] += dt*s->atoms->p[iOff][2]*invMass;
		}
	}
	#ifdef __pointerchain
	#pragma pointerchain region end
	#endif
}

/// Calculates total kinetic and potential energy across all tasks.  The
/// local potential energy is a by-product of the force routine.
void kineticEnergy(SimFlat* s, int ntoms_present)
{
	real_t eLocal[2];
	real_t kenergy = 0.0;
	eLocal[0] = s->ePotential;
	eLocal[1] = 0;

	#ifdef __pointerchain
	#pragma pointerchain declare(s->species{SpeciesData*}, s->boxes->nTotalBoxes{int}, s->boxes->nbrBoxes{int**}, s->boxes->nLocalBoxes{int})
	#pragma pointerchain declare(s->atoms->r{real3*}, s->atoms->p{real3*}, s->atoms->f{real3*}, s->atoms->U{real_t*})
	#pragma pointerchain declare(s->boxes->nAtoms{int*}, s->atoms->iSpecies{int*}, s->atoms->gid{int*})
	#endif

	int fSize = s->boxes->nTotalBoxes*MAXATOMS;

	#ifdef __pointerchain
	#pragma pointerchain region begin
	#endif
	#pragma acc parallel loop gang independent present(s->boxes->nAtoms[0:s->boxes->nTotalBoxes], \
		s->atoms->p[0:fSize][0:3], s->atoms->iSpecies[0:fSize], s->species[0:1]) \
		if(ntoms_present)
	for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
	{
		#pragma acc loop vector independent
		for (int ii=0; ii<s->boxes->nAtoms[iBox]; ii++)
		{
			int iOff = MAXATOMS * iBox + ii;
			int iSpecies = s->atoms->iSpecies[iOff];
			real_t invMass = 0.5/s->species[iSpecies].mass;
			kenergy += ( s->atoms->p[iOff][0] * s->atoms->p[iOff][0] +
						s->atoms->p[iOff][1] * s->atoms->p[iOff][1] +
						s->atoms->p[iOff][2] * s->atoms->p[iOff][2] )*invMass;
		}
	}
	#ifdef __pointerchain
	#pragma pointerchain region end
	#endif

	eLocal[1] = kenergy;

	real_t eSum[2];
	startTimer(commReduceTimer);
	addRealParallel(eLocal, eSum, 2);
	stopTimer(commReduceTimer);

	s->ePotential = eSum[0];
	s->eKinetic = eSum[1];
}

/// \details
/// This function provides one-stop shopping for the sequence of events
/// that must occur for a proper exchange of halo atoms after the atom
/// positions have been updated by the integrator.
///
/// - updateLinkCells: Since atoms have moved, some may be in the wrong
///   link cells.
/// - haloExchange (atom version): Sends atom data to remote tasks. 
/// - sort: Sort the atoms.
///
/// \see updateLinkCells
/// \see initAtomHaloExchange
/// \see sortAtomsInCell
void redistributeAtomsInInit(SimFlat* sim) {
//	startTimer(updateLinkCellsTimer);
	updateLinkCells(sim->boxes, sim->atoms);
//	stopTimer(updateLinkCellsTimer);

	startTimer(atomHaloTimer);
	haloExchange(sim->atomExchange, sim);
	stopTimer(atomHaloTimer);

//	startTimer(sortAtomsTimer);
#pragma omp parallel for
	for (int ii = 0; ii < sim->boxes->nTotalBoxes; ++ii)
		sortAtomsInCell(sim->atoms, sim->boxes, ii);
//	stopTimer(sortAtomsTimer);

}

void redistributeAtoms(SimFlat* sim)
{
	updateHost(sim);

//	startTimer(updateLinkCellsTimer);
	updateLinkCells(sim->boxes, sim->atoms);
//	stopTimer(updateLinkCellsTimer);

	startTimer(atomHaloTimer);
	haloExchange(sim->atomExchange, sim);
	stopTimer(atomHaloTimer);

//	startTimer(sortAtomsTimer);
	#pragma omp parallel for
	for (int ii=0; ii<sim->boxes->nTotalBoxes; ++ii)
		sortAtomsInCell(sim->atoms, sim->boxes, ii);
//	stopTimer(sortAtomsTimer);

	updateGPU(sim);
}
