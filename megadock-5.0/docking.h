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
//  Class Name : Docking
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Docking_h
#define Docking_h 1

#include "parallel.h"
#include "receptor.h"
#include "ligand.h"

using namespace std;

template<class P, class F> class Docking
{
private:
    Docking(Docking &c) {}
    const Docking & operator=(const Docking &c);
protected:
    Parallel  *_parallel;
    P *_parameter;
    Receptor<P>  *_receptor;
    Ligand<P>    *_ligand;
    F    *_fft_process;
    int       _Num_grid;
    float     *_Grid_coord;
    float     *_Fwork;
    int       *_Iwork;
    size_t       _Memfw;
    size_t       _Memiw;
    virtual void  maxsize_voxel() = 0;
    virtual void  alloc_array(const int &maxatom, const int &nag, const size_t &ng3) = 0;
public:
    Docking(Parallel *pparallel,P *pparameter,
            Receptor<P> *rreceptor,Ligand<P> *rligand)
        : _parallel(pparallel),_parameter(pparameter),
          _receptor(rreceptor),_ligand(rligand) {
#ifdef DEBUG
        cout << "Constructing Docking.\n";
#endif
    }
    virtual ~Docking() {
#ifdef DEBUG
        cout << "Destructing Docking.\n";
#endif
        delete [] _Grid_coord;
        delete [] _Fwork;
        delete [] _Iwork;
        delete _fft_process;
    }
    virtual void  initialize() = 0;
    virtual void  rec_init() = 0;
    virtual void  dockz() = 0;
    virtual void  dock_memory_free() = 0;
    virtual void  output() = 0;
    virtual void  output_detail() = 0; // for analysis
    virtual void  output_calc_time_log() = 0; // for analysis
};

#endif
