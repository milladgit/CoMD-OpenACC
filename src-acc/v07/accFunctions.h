/*
 * accFunctions.h
 *
 *  Created on: Dec 17, 2016
 *      Author: millad
 */

#ifndef ACCFUNCTIONS_H_
#define ACCFUNCTIONS_H_

#include "CoMDTypes.h"
#include "linkCells.h"
#include "parallel.h"

#include <openacc.h>

static inline
void SetupGpu(int deviceId)
{
	acc_set_device_num(deviceId, acc_device_nvidia);

	char hostname[256];
	gethostname(hostname, sizeof(hostname));

	printf("Host %s using NVIDIA GPU %i\n\n", hostname, deviceId);
}


static inline 
void transferToGPU(SimFlat *s) {
	int nNbrBoxes = 27;

	#ifdef __pointerchain
	#pragma pointerchain declare(s->species{SpeciesData*}, s->boxes->nTotalBoxes{int}, s->boxes->nbrBoxes{int**})
	#pragma pointerchain declare(s->atoms->rx{real_t*}, s->atoms->ry{real_t*}, s->atoms->rz{real_t*})
	#pragma pointerchain declare(s->atoms->px{real_t*}, s->atoms->py{real_t*}, s->atoms->pz{real_t*})
	#pragma pointerchain declare(s->atoms->fx{real_t*}, s->atoms->fy{real_t*}, s->atoms->fz{real_t*})
	#pragma pointerchain declare(s->atoms->U{real_t*})
	#pragma pointerchain declare(s->boxes->nAtoms{int*}, s->atoms->iSpecies{int*}, s->atoms->gid{int*})
	#endif


	/*
	 * Converting s->boxes->nbrBoxes[i][j] to a flat array
	 */
//	int* s_boxes_nbrBoxes = malloc(sizeof(int) * s->boxes->nTotalBoxes * nNbrBoxes);
//	s->boxes->nbrBoxes_array = s_boxes_nbrBoxes;
//	for(int i=0;i<s->boxes->nTotalBoxes;i++)
//		memcpy(&s_boxes_nbrBoxes[i * nNbrBoxes], s->boxes->nbrBoxes[i], nNbrBoxes * sizeof(int));
//	#pragma acc enter data copyin(s_boxes_nbrBoxes[0:s->boxes->nTotalBoxes * nNbrBoxes])


	int maxTotalAtoms = MAXATOMS * s->boxes->nTotalBoxes;
	int rf_size = maxTotalAtoms * sizeof(real3);
	int U_size = maxTotalAtoms * sizeof(real_t);

	#ifdef __pointerchain
	#pragma pointerchain region begin
	#endif
	#pragma acc enter data copyin(s->atoms->rx[0:maxTotalAtoms], s->atoms->ry[0:maxTotalAtoms], s->atoms->rz[0:maxTotalAtoms])
	#pragma acc enter data copyin(s->atoms->px[0:maxTotalAtoms], s->atoms->py[0:maxTotalAtoms], s->atoms->pz[0:maxTotalAtoms])
	#pragma acc enter data copyin(s->atoms->fx[0:maxTotalAtoms], s->atoms->fy[0:maxTotalAtoms], s->atoms->fz[0:maxTotalAtoms])
	#pragma acc enter data copyin(s->atoms->U[0:maxTotalAtoms], s->boxes->nAtoms[0:s->boxes->nTotalBoxes], s->atoms->iSpecies[0:maxTotalAtoms], s->atoms->gid[0:maxTotalAtoms])
	#pragma acc enter data copyin(s->species[0:1], s->boxes->nbrBoxes[0:s->boxes->nTotalBoxes][0:nNbrBoxes])
	#ifdef __pointerchain
	#pragma pointerchain region end
	#endif
}



static inline 
void updateHost(SimFlat *s) {

	int maxTotalAtoms = MAXATOMS * s->boxes->nTotalBoxes;
	int rf_size = maxTotalAtoms * sizeof(real3);
	int U_size = maxTotalAtoms * sizeof(real_t);

	#ifdef __pointerchain
	#pragma pointerchain declare(s->boxes->nTotalBoxes{int})
	#pragma pointerchain declare(s->atoms->rx{real_t*}, s->atoms->ry{real_t*}, s->atoms->rz{real_t*})
	#pragma pointerchain declare(s->atoms->px{real_t*}, s->atoms->py{real_t*}, s->atoms->pz{real_t*})
	#pragma pointerchain declare(s->atoms->fx{real_t*}, s->atoms->fy{real_t*}, s->atoms->fz{real_t*})
	#pragma pointerchain declare(s->atoms->U{real_t*})
	#pragma pointerchain declare(s->boxes->nAtoms{int*}, s->atoms->iSpecies{int*}, s->atoms->gid{int*})
	#pragma pointerchain declare(s->boxes->localMin{real_t*}, s->boxes->localMax{real_t*}, s->boxes->invBoxSize{real_t*}, s->boxes->gridSize{int*})
	#endif


	#ifdef __pointerchain
	#pragma pointerchain region begin
	#endif
	#pragma acc update self(s->atoms->rx[0:maxTotalAtoms], s->atoms->ry[0:maxTotalAtoms], s->atoms->rz[0:maxTotalAtoms])
	#pragma acc update self(s->atoms->px[0:maxTotalAtoms], s->atoms->py[0:maxTotalAtoms], s->atoms->pz[0:maxTotalAtoms])
	#pragma acc update self(s->atoms->fx[0:maxTotalAtoms], s->atoms->fy[0:maxTotalAtoms], s->atoms->fz[0:maxTotalAtoms])
	#pragma acc update self(s->atoms->U[0:maxTotalAtoms], s->boxes->nAtoms[0:s->boxes->nTotalBoxes], s->atoms->gid[0:maxTotalAtoms])
//	#pragma acc update self (s->boxes->localMin[0:3], s->boxes->localMax[0:3], s->boxes->invBoxSize[0:3], s->boxes->gridSize[0:3])
	#ifdef __pointerchain
	#pragma pointerchain region end
	#endif
}

static inline 
void updateGPU(SimFlat *s) {

	int maxTotalAtoms = MAXATOMS * s->boxes->nTotalBoxes;
	int rf_size = maxTotalAtoms * sizeof(real3);
	int U_size = maxTotalAtoms * sizeof(real_t);

	#ifdef __pointerchain
	#pragma pointerchain declare(s->boxes->nTotalBoxes{int})
	#pragma pointerchain declare(s->atoms->rx{real_t*}, s->atoms->ry{real_t*}, s->atoms->rz{real_t*})
	#pragma pointerchain declare(s->atoms->px{real_t*}, s->atoms->py{real_t*}, s->atoms->pz{real_t*})
	#pragma pointerchain declare(s->atoms->fx{real_t*}, s->atoms->fy{real_t*}, s->atoms->fz{real_t*})
	#pragma pointerchain declare(s->atoms->U{real_t*})
	#pragma pointerchain declare(s->boxes->nAtoms{int*}, s->atoms->iSpecies{int*}, s->atoms->gid{int*})
	#pragma pointerchain declare(s->boxes->localMin{real_t*}, s->boxes->localMax{real_t*}, s->boxes->invBoxSize{real_t*}, s->boxes->gridSize{int*})
	#endif

	#ifdef __pointerchain
	#pragma pointerchain region begin
	#endif
	#pragma acc update device(s->atoms->rx[0:maxTotalAtoms], s->atoms->ry[0:maxTotalAtoms], s->atoms->rz[0:maxTotalAtoms])
	#pragma acc update device(s->atoms->px[0:maxTotalAtoms], s->atoms->py[0:maxTotalAtoms], s->atoms->pz[0:maxTotalAtoms])
	#pragma acc update device(s->atoms->fx[0:maxTotalAtoms], s->atoms->fy[0:maxTotalAtoms], s->atoms->fz[0:maxTotalAtoms])
	#pragma acc update device(s->atoms->U[0:maxTotalAtoms], s->boxes->nAtoms[0:s->boxes->nTotalBoxes], s->atoms->gid[0:maxTotalAtoms])
//	#pragma acc update device (s->boxes->localMin[0:3], s->boxes->localMax[0:3], s->boxes->invBoxSize[0:3], s->boxes->gridSize[0:3])
	#ifdef __pointerchain
	#pragma pointerchain region end
	#endif
}


#endif /* ACCFUNCTIONS_H_ */
