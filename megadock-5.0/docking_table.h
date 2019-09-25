/*
 * Copyright (C) 2014 Tokyo Institute of Technology
 *
 *
 * This file is part of MEGADOCK.
 * MEGADOCK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MEGADOCK is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MEGADOCK.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : DockingTable
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef DockingTable_h
#define DockingTable_h 1

#include "docking.h"
#include "exec_logger.h"
#include "parameter_table.h"
#include "fft_process_table.h"

using namespace std;

class DockingTable : public Docking<ParameterTable, FFTProcessTable>
{
private:
    ExecLogger   *_exec_logger;
    float     *_Mol_coord;
protected:
    virtual void  maxsize_voxel();
    virtual void  alloc_array(const int &maxatom, const int &nag, const size_t &ng3);
    virtual void  create_voxel(Protein<ParameterTable> *rprotein);
    virtual void  ligand_rotationz(float *theta);
public:
    DockingTable(ExecLogger *pexec_logger,Parallel *pparallel,ParameterTable *pparameter,
            Receptor<ParameterTable> *rreceptor,Ligand<ParameterTable> *rligand)
        : _exec_logger(pexec_logger),Docking(pparallel,pparameter,rreceptor,rligand) {
#ifdef DEBUG
        cout << "Constructing DockingTable.\n";
#endif
    }
    virtual ~DockingTable() {
#ifdef DEBUG
        cout << "Destructing DockingTable.\n";
#endif
    }
    virtual void  initialize();
    virtual void  rec_init();
    virtual void  dockz();
    virtual void  dock_memory_free();
    virtual void  output();
    virtual void  output_detail(); // for analysis
    virtual void  output_calc_time_log(); // for analysis
};

#endif
