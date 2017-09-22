/*
    This file is part of the LibTri library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Paola Magillo (paola.magillo@unige.it)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The LibTri library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The LibTri library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the LibTri library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SORTING_H
#define SORTING_H

typedef struct {
    int v1,v2;
    int t;
} aux;

static int cmp_aux( const void *p, const void *q )
{
  aux *a = (aux *) p ;
  aux *b = (aux *) q ;
  if ( a->v1 > b->v1 ) return 1 ;
  if ( a->v1 < b->v1 ) return -1 ;
  if ( a->v2 > b->v2 ) return 1;
  if ( a->v2 < b->v2 ) return -1;
  return 0;
}

#endif // SORTING_H
