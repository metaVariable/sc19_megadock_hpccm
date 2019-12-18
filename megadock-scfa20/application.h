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
//  Class Name : Application
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Application_h
#define Application_h 1

#include "cpu_time.h"
#include "exec_logger.h"
#include "control_pdb.h"
#include "control_table.h"

class Application
{
private:
    Application(Application &c) {}
    const Application & operator=(const Application &c);

    int nproc2;
    int device_count_gpu;
    Parallel  **_parallels;
    ExecLogger   **_exec_loggers;
    ControlTable   **_controls;
    ParameterTable **_parameters;

public:
    Application(const int nproc2) : nproc2(nproc2) {}
    virtual ~Application() {
#pragma omp parallel for
        for (int i = 0; i < nproc2; i++) {
            delete _exec_loggers[i];
            delete _controls[i];
            delete _parallels[i];
            delete _parameters[i];
        }

        delete [] _parallels;
        delete [] _exec_loggers;
        delete [] _controls;
        delete [] _parameters;
    }
    virtual void initialize();
    virtual int application(int argc, char *argv[], int myid2);
};

#endif
