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
	#pragma pointerchain declare(s->atoms->px{real_t*}, s->atoms->py{real_t*}, s->atoms->pz{real_t*})
	#pragma pointerchain declare(s->atoms->fx{real_t*}, s->atoms->fy{real_t*}, s->atoms->fz{real_t*})
	#pragma pointerchain declare(s->boxes->nAtoms{int*}, s->atoms->iSpecies{int*}, s->atoms->gid{int*})
	#endif

	int fSize = s->boxes->nTotalBoxes*MAXATOMS;

	#ifdef __pointerchain
	#pragma pointerchain region begin
	#endif
	#pragma acc parallel loop num_gangs(s->gangCount) vector_length(s->vectorCount) \
		gang vector collapse(2) independent \
		present(s->boxes->nAtoms[0:nBoxes], \
				s->atoms->px[0:fSize], s->atoms->py[0:fSize], s->atoms->pz[0:fSize], \
				s->atoms->fx[0:fSize], s->atoms->fy[0:fSize], s->atoms->fz[0:fSize])
	for (int iBox=0; iBox<nBoxes; iBox++)
	{
		for (int ii=0; ii<MAXATOMS; ii++)
		{
			if(ii >= s->boxes->nAtoms[iBox])
				continue;

			int iOff = MAXATOMS * iBox + ii;

			s->atoms->px[iOff] += dt * s->atoms->fx[iOff];
			s->atoms->py[iOff] += dt * s->atoms->fy[iOff];
			s->atoms->pz[iOff] += dt * s->atoms->fz[iOff];
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
	#pragma pointerchain declare(s->atoms->rx{real_t*}, s->atoms->ry{real_t*}, s->atoms->rz{real_t*})
	#pragma pointerchain declare(s->atoms->px{real_t*}, s->atoms->py{real_t*}, s->atoms->pz{real_t*})
	#pragma pointerchain declare(s->boxes->nAtoms{int*}, s->atoms->iSpecies{int*}, s->atoms->gid{int*})
	#endif

	int fSize = s->boxes->nTotalBoxes*MAXATOMS;

	#ifdef __pointerchain
	#pragma pointerchain region begin
	#endif
	#pragma acc parallel loop num_gangs(s->gangCount) vector_length(s->vectorCount) \
		gang vector collapse(2) independent \
		present(s->boxes->nAtoms[0:nBoxes], s->atoms->iSpecies[0:fSize], s->species[0:1], \
				s->atoms->px[0:fSize], s->atoms->py[0:fSize], s->atoms->pz[0:fSize], \
				s->atoms->rx[0:fSize], s->atoms->ry[0:fSize], s->atoms->rz[0:fSize])
	for (int iBox=0; iBox<nBoxes; iBox++)
	{
		for (int ii=0; ii<MAXATOMS; ii++)
		{
			if(ii >= s->boxes->nAtoms[iBox])
				continue;

			int iOff = MAXATOMS * iBox + ii;
			int iSpecies = s->atoms->iSpecies[iOff];
			real_t invMass = 1.0/s->species[iSpecies].mass;
			s->atoms->rx[iOff] += dt * s->atoms->px[iOff] * invMass;
			s->atoms->ry[iOff] += dt * s->atoms->py[iOff] * invMass;
			s->atoms->rz[iOff] += dt * s->atoms->pz[iOff] * invMass;
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
	real_t penergy = 0.0;
	eLocal[0] = 0;
	eLocal[1] = 0;

	#ifdef __pointerchain
	#pragma pointerchain declare(s->species{SpeciesData*}, s->boxes->nTotalBoxes{int}, s->boxes->nbrBoxes{int**}, s->boxes->nLocalBoxes{int})
	#pragma pointerchain declare(s->atoms->rx{real_t*}, s->atoms->ry{real_t*}, s->atoms->rz{real_t*})
	#pragma pointerchain declare(s->atoms->px{real_t*}, s->atoms->py{real_t*}, s->atoms->pz{real_t*})
	#pragma pointerchain declare(s->atoms->U{real_t*})
	#pragma pointerchain declare(s->boxes->nAtoms{int*}, s->atoms->iSpecies{int*}, s->atoms->gid{int*})
	#endif

	int fSize = s->boxes->nTotalBoxes*MAXATOMS;

	#ifdef __pointerchain
	#pragma pointerchain region begin
	#endif
	#pragma acc parallel loop \
		gang vector reduction(+:kenergy) collapse(2) independent \
		present(s->boxes->nAtoms[0:s->boxes->nTotalBoxes], \
				s->atoms->px[0:fSize], s->atoms->py[0:fSize], s->atoms->pz[0:fSize],\
				s->atoms->iSpecies[0:fSize], s->species[0:1]) \
		if(ntoms_present)
	for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
	{
		for (int ii=0; ii<MAXATOMS; ii++)
		{
			if(ii >= s->boxes->nAtoms[iBox])
				continue;

			int iOff = MAXATOMS * iBox + ii;
			int iSpecies = s->atoms->iSpecies[iOff];
			real_t invMass = 0.5/s->species[iSpecies].mass;
			real_t px = s->atoms->px[iOff];
			real_t py = s->atoms->py[iOff];
			real_t pz = s->atoms->pz[iOff];
			kenergy += (px*px + py*py + pz*pz) * invMass;
		}
	}

	#pragma acc parallel loop \
		gang vector reduction(+:penergy) collapse(2) independent \
		present(s->boxes->nAtoms[0:s->boxes->nTotalBoxes], s->atoms->U[0:fSize]) \
		if(ntoms_present)
	for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
	{
		for (int ii=0; ii<MAXATOMS; ii++)
		{
			if(ii >= s->boxes->nAtoms[iBox])
				continue;

			int iOff = MAXATOMS * iBox + ii;
			penergy += s->atoms->U[iOff];
		}
	}

	#ifdef __pointerchain
	#pragma pointerchain region end
	#endif

	eLocal[0] = penergy;
	eLocal[1] = kenergy;

	s->ePotential = penergy;

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
