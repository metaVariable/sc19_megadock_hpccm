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
//  Class Name : Control
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Control_h
#define Control_h 1

#include "parallel.h"
#include "receptor.h"
#include "ligand.h"
#include "docking_pdb.h"
#include "docking_table.h"

using namespace std;

template<class P, class D> class Control
{
private:
    Control(Control &c) {}
    const Control & operator=(const Control &c);
protected:
    Parallel  *_parallel;
    P *_parameter;
    Receptor<P>  *_receptor;
    Ligand<P>    *_ligand;
    D            *_docking;
    virtual void  gridtable_11base_normal(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_13base_normal(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_07base_normal(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_fftw_custom(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_cufft_custom(int &ngrid,vector<int> &ngrid_table);
    virtual void  autogridr(const int &ngrid,vector<int> &ngrid_table) = 0;
    virtual void  autogridl(const int &ngrid,vector<int> &ngrid_table) = 0;
    virtual void  checkgridr() = 0;
    virtual void  checkgridl() = 0;
public:
    Control(Parallel *pparallel) : _parallel(pparallel) {
#ifdef DEBUG
        cout << "Constructing Control.\n";
#endif
    }
    Control(Parallel *pparallel, P *pparameter) : _parallel(pparallel), _parameter(pparameter) {
#ifdef DEBUG
        cout << "Constructing Control.\n";
#endif
    }
    virtual ~Control() {
#ifdef DEBUG
        cout << "Destructing Control.\n";
#endif
    }
    virtual void  execute() = 0;
};

#endif
