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
//  Class Name : ParameterPDB
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef ParameterPDB_h
#define ParameterPDB_h 1

#include "parameter.h"

using namespace std;

class ParameterPDB : public Parameter
{
private:
    friend class          ControlPDB;
    friend class          DockingPDB;
    friend class          FFTProcessPDB;
    string            _RecPDB_file;
    string            _LigPDB_file;
    int               _IO_flag[3];

protected:
    virtual void          pdb_step();

public:
    ParameterPDB(Parallel *pparallel) : Parameter(pparallel) {
#ifdef DEBUG
        cout << "Constructing ParameterPDB.\n";
#endif
    }
    virtual           ~ParameterPDB() {
#ifdef DEBUG
        cout << "Destructing ParameterPDB.\n";
#endif
    }
    virtual void          process_args(int argc,char *argv[]);
};

#endif
