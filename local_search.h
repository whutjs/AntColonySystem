/*

       AAAA    CCCC   OOOO   TTTTTT   SSSSS  PPPPP
      AA  AA  CC     OO  OO    TT    SS      PP  PP
      AAAAAA  CC     OO  OO    TT     SSSS   PPPPP
      AA  AA  CC     OO  OO    TT        SS  PP
      AA  AA   CCCC   OOOO     TT    SSSSS   PP

######################################################
##########    ACO algorithms for the TSP    ##########
######################################################

      Version: 1.0
      File:    ls.h
      Author:  Thomas Stuetzle
      Purpose: header file for local search routines
      Check:   README and gpl.txt
      Copyright (C) 1999  Thomas Stuetzle
*/

/***************************************************************************

    Program's name: acotsp

    Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the 
    symmetric TSP 

    Copyright (C) 2004  Thomas Stuetzle

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    email: stuetzle no@spam ulb.ac.be
    mail address: Universite libre de Bruxelles
                  IRIDIA, CP 194/6
                  Av. F. Roosevelt 50
                  B-1050 Brussels
		  Belgium

***************************************************************************/

#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

#include "ants.h"
#include <math.h>

extern int ls_flag;

extern int nn_ls; 

extern bool dlb_flag; 

extern unsigned short *pos;               /* positions of cities in tour */ 
extern bool *dlb;               /* vector containing don't look bits */ 
extern unsigned short *h_tour;            /* help vector for performing exchange move */ 
extern unsigned short *hh_tour;           /* help vector for performing exchange move */ 
extern unsigned short *random_vector;
/*
	n=nodes.size();
*/
//void three_opt_first(std::vector<unsigned short>& tour, const unsigned short& n, 
//					 const bool* visited,
//					 const unsigned short& ghost_cnt);

void three_opt_first(std::vector<unsigned short>& tour, const unsigned short& n);

void init_ls();
void free_ls();

#endif