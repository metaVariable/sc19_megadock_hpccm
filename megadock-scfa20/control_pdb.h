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
//  Class Name : ControlPDB
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef ControlPDB_h
#define ControlPDB_h 1

#include "control.h"
#include "cpu_time.h"
#include "parameter_pdb.h"
#include "docking.h"

using namespace std;

class ControlPDB : public Control<ParameterPDB, DockingPDB>
{
private:
    CPUTime   *_cputime;
protected:
    virtual void  autogridr(const int &ngrid,vector<int> &ngrid_table);
    virtual void  autogridl(const int &ngrid,vector<int> &ngrid_table);
    virtual void  checkgridr();
    virtual void  checkgridl();
public:
    ControlPDB(CPUTime *pcputime,Parallel *pparallel)
        : _cputime(pcputime),Control<ParameterPDB, DockingPDB>(pparallel) {
#ifdef DEBUG
        cout << "Constructing ControlPDB.\n";
#endif
    }
    virtual ~ControlPDB() {
#ifdef DEBUG
        cout << "Destructing ControlPDB.\n";
#endif
    }
    virtual void  initialize(int argc,char *argv[]);
    virtual void  execute();
};

#endif
