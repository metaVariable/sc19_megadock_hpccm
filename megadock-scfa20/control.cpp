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

#include "control.h"

//============================================================================//
template<class P, class D> void Control<P, D>::gridtable_11base_normal(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = {2,3,4,5,6,7,8,9,10,11,12,14,15,16,18,20,21,22,24,25,27,
                        28,30,32,33,35,36,40,42,44,45,48,49,50,54,55,56,60,63,
                        64,66,70,72,75,77,80,81,84,88,90,96,98,99,100,105,108,
                        110,112,120,121,125,126,128,132,135,140,144,147,150,
                        154,160,162,165,168,175,176,180,189,192,196,198,200,
                        210,216,220,224,225,231,240,242,243,245,250,252,256,
                        264,270,275,280,288,294,297,300,308,315,320,324,330,
                        336,343,350,352,360,363,375,378,384,385,392,396,400,
                        405,420,432,440,441,448,450,462,480,484,486,490,495,
                        500,504,512,525,528,539,540,550,560,567,576,588,594,
                        600,605,616,625,630,640,648,660,672,675,686,693,700,
                        704,720,726,729,735,750,756,768,770,784,792,800,810,
                        825,840,847,864,875,880,882,891,896,900,924,945,960,
                        968,972,980,990,1000,1008,1024
                       };

    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}

//============================================================================//
template<class P, class D> void Control<P, D>::gridtable_13base_normal(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = { 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 ,
                         12 , 13 , 14 , 15 , 16 , 18 , 20 , 21 , 22 , 24 ,
                         25 , 26 , 27 , 28 , 30 , 32 , 33 , 35 , 36 , 39 ,
                         40 , 42 , 44 , 45 , 48 , 49 , 50 , 52 , 54 , 55 ,
                         56 , 60 , 63 , 64 , 65 , 66 , 70 , 72 , 75 , 77 ,
                         78 , 80 , 81 , 84 , 88 , 90 , 91 , 96 , 98 , 99 ,
                         100 , 104 , 105 , 108 , 110 , 112 , 117 , 120 , 121 , 125 ,
                         126 , 128 , 130 , 132 , 135 , 140 , 143 , 144 , 147 , 150 ,
                         154 , 156 , 160 , 162 , 165 , 168 , 169 , 175 , 176 , 180 ,
                         182 , 189 , 192 , 195 , 196 , 198 , 200 , 208 , 210 , 216 ,
                         220 , 224 , 225 , 231 , 234 , 240 , 242 , 243 , 245 , 250 ,
                         252 , 256 , 260 , 264 , 270 , 273 , 275 , 280 , 286 , 288 ,
                         294 , 297 , 300 , 308 , 312 , 315 , 320 , 324 , 325 , 330 ,
                         336 , 338 , 343 , 350 , 351 , 352 , 360 , 363 , 364 , 375 ,
                         378 , 384 , 385 , 390 , 392 , 396 , 400 , 405 , 416 , 420 ,
                         429 , 432 , 440 , 441 , 448 , 450 , 455 , 462 , 468 , 480 ,
                         484 , 486 , 490 , 495 , 500 , 504 , 507 , 512 , 520 , 525 ,
                         528 , 539 , 540 , 546 , 550 , 560 , 567 , 572 , 576 , 585 ,
                         588 , 594 , 600 , 605 , 616 , 624 , 625 , 630 , 637 , 640 ,
                         648 , 650 , 660 , 672 , 675 , 676 , 686 , 693 , 700 , 702 ,
                         704 , 715 , 720 , 726 , 728 , 729 , 735 , 750 , 756 , 768 ,
                         770 , 780 , 784 , 792 , 800 , 810 , 819 , 825 , 832 , 840 ,
                         845 , 847 , 858 , 864 , 875 , 880 , 882 , 891 , 896 , 900 ,
                         910 , 924 , 936 , 945 , 960 , 968 , 972 , 975 , 980 , 990 ,
                         1000 , 1001 , 1008 , 1014 , 1024
                       };

    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}

//============================================================================//
template<class P, class D> void Control<P, D>::gridtable_07base_normal(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = { 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 12 ,
                         14 , 15 , 16 , 18 , 20 , 21 , 24 , 25 , 27 , 28 ,
                         30 , 32 , 35 , 36 , 40 , 42 , 45 , 48 , 49 , 50 ,
                         54 , 56 , 60 , 63 , 64 , 70 , 72 , 75 , 80 , 81 ,
                         84 , 90 , 96 , 98 , 100 , 105 , 108 , 112 , 120 , 125 ,
                         126 , 128 , 135 , 140 , 144 , 147 , 150 , 160 , 162 , 168 ,
                         175 , 180 , 189 , 192 , 196 , 200 , 210 , 216 , 224 , 225 ,
                         240 , 243 , 245 , 250 , 252 , 256 , 270 , 280 , 288 , 294 ,
                         300 , 315 , 320 , 324 , 336 , 343 , 350 , 360 , 375 , 378 ,
                         384 , 392 , 400 , 405 , 420 , 432 , 441 , 448 , 450 , 480 ,
                         486 , 490 , 500 , 504 , 512 , 525 , 540 , 560 , 567 , 576 ,
                         588 , 600 , 625 , 630 , 640 , 648 , 672 , 675 , 686 , 700 ,
                         720 , 729 , 735 , 750 , 756 , 768 , 784 , 800 , 810 , 840 ,
                         864 , 875 , 882 , 896 , 900 , 945 , 960 , 972 , 980 , 1000 ,
                         1008 , 1024
                       };


    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}

//============================================================================//
template<class P, class D> void Control<P, D>::gridtable_fftw_custom(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = { 8,10,16,18,21,22,24,25,32,33,35,
                         36,39,40,42,44,45,49,50,52,55,56,60,64,70,75,78,84,91,
                         100,105,110,117,126,130,140,147,150,154,156,162,165,169,175,189,200,210,216,220,225,
                         231 , 234 , 240 , 242 , 243 , 245 , 250 ,
                         252 , 256 , 260 , 264 , 270 , 273 , 275 , 280 , 286 , 288 ,
                         294 , 297 , 300 , 308 , 312 , 315 , 320 , 324 , 325 , 330 ,
                         336 , 338 , 343 , 350 , 351 , 352 , 360 , 363 , 364 , 375 ,
                         378 , 384 , 385 , 390 , 392 , 396 , 400 , 405 , 416 , 420 ,
                         429 , 432 , 440 , 441 , 448 , 450 , 455 , 462 , 468 , 480 ,
                         484 , 486 , 490 , 495 , 500 , 504 , 507 , 512 , 520 , 525 ,
                         528 , 539 , 540 , 546 , 550 , 560 , 567 , 572 , 576 , 585 ,
                         588 , 594 , 600 , 605 , 616 , 624 , 625 , 630 , 637 , 640 ,
                         648 , 650 , 660 , 672 , 675 , 676 , 686 , 693 , 700 , 702 ,
                         704 , 715 , 720 , 726 , 728 , 729 , 735 , 750 , 756 , 768 ,
                         770 , 780 , 784 , 792 , 800 , 810 , 819 , 825 , 832 , 840 ,
                         845 , 847 , 858 , 864 , 875 , 880 , 882 , 891 , 896 , 900 ,
                         910 , 924 , 936 , 945 , 960 , 968 , 972 , 975 , 980 , 990 ,
                         1000 , 1001 , 1008 , 1014 , 1024
                       }; //selected the just values to have faster FFT calculation time (in n < 226)


    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}

//============================================================================//
template<class P, class D> void Control<P, D>::gridtable_cufft_custom(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = { 16,32,36,40,42,48,50,64,72,80,81,84,96,98,100,108,112,128,144,160,162,168,200,216,224,240 ,
                         243 , 245 , 250 , 252 , 256 , 270 , 280 , 288 , 294 ,
                         300 , 315 , 320 , 324 , 336 , 343 , 350 , 360 , 375 , 378 ,
                         384 , 392 , 400 , 405 , 420 , 432 , 441 , 448 , 450 , 480 ,
                         486 , 490 , 500 , 504 , 512 , 525 , 540 , 560 , 567 , 576 ,
                         588 , 600 , 625 , 630 , 640 , 648 , 672 , 675 , 686 , 700 ,
                         720 , 729 , 735 , 750 , 756 , 768 , 784 , 800 , 810 , 840 ,
                         864 , 875 , 882 , 896 , 900 , 945 , 960 , 972 , 980 , 1000 ,
                         1008 , 1024
                       }; //selected the just values to have faster FFT calculation time (in n < 226)


    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}

// explicit instantiations
template class Control<ParameterPDB, DockingPDB>;
template class Control<ParameterTable, DockingTable>;
