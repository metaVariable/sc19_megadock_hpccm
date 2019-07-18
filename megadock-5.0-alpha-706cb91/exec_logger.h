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
//  Class Name : ExecLogger
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Exec_logger_h
#define Exec_logger_h 1

#include <stdio.h>
#include <string>
#include <iostream>

#include "cpu_time.h"

using namespace std;

class ExecLogger
{
private:
    ExecLogger(ExecLogger &c) {}
    const ExecLogger & operator=(const ExecLogger &c);

public:
    CPUTime *_cputime;
    long double mb;

#ifdef CUFFT
    size_t devmem_free, devmem_total, devmem_use;
#endif

    string _RLOut_file;
    int _Num_fft_flag;
    string rec_filename, lig_filename;
    float rec_max_size, lig_max_size;
    float rec_voxel_size, lig_voxel_size;
    int rec_num_grid, lig_num_grid;
    float grid_width;

    ExecLogger() {
#ifdef DEBUG
        cout << "Constructing ExecLogger.\n";
#endif
        _cputime = new CPUTime();
    }
    virtual   ~ExecLogger()      {
#ifdef DEBUG
        cout << "Destructing ExecLogger\n";
#endif
        delete _cputime;
    }
    virtual void      initialize();
    virtual void      output(const int myid2);
    virtual void      record_malloc(const int &size) {
        _cputime->record_malloc(size);
    }
    virtual void      record_free(const int &size) {
        _cputime->record_free(size);
    }
};

#endif
