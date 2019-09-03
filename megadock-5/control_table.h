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
//  Class Name : ControlTable
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef ControlTable_h
#define ControlTable_h 1

#include "control.h"
#include "exec_logger.h"
#include "parameter_table.h"
#include "docking_table.h"

using namespace std;

class ControlTable : public Control<ParameterTable, DockingTable>
{
private:
    ExecLogger   *_exec_logger;
protected:
    virtual void  autogridr(const int &ngrid,vector<int> &ngrid_table);
    virtual void  autogridl(const int &ngrid,vector<int> &ngrid_table);
    virtual void  checkgridr();
    virtual void  checkgridl();
public:
    ControlTable(ExecLogger *pexec_logger,Parallel *pparallel,ParameterTable *pparameter)
        : _exec_logger(pexec_logger),Control<ParameterTable, DockingTable>(pparallel,pparameter) {
#ifdef DEBUG
        cout << "Constructing ControlTable.\n";
#endif
    }
    virtual ~ControlTable() {
#ifdef DEBUG
        cout << "Destructing ControlTable.\n";
#endif
    }
    virtual void  initialize(bool verbose);
    virtual void  prepare(string rec_file, string lig_file, string out_file);
    virtual void  execute();
    virtual string input_file() {
        return _parameter->_Table_file;
    }
};

#endif
