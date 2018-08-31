/*
 * Copyright (C) 1999-2006 Matthias Seeger
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 */
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Definition class Interval
 * ------------------------------------------------------------------- */

#include "lhotse/global.h"
#include "lhotse/Interval.h"

template<> bool DefIVal<char>::isInit(false);
template<> bool DefIVal<uchar>::isInit(false);
template<> bool DefIVal<int>::isInit(false);
template<> bool DefIVal<uint>::isInit(false);
template<> bool DefIVal<long>::isInit(false);
template<> bool DefIVal<unsigned long>::isInit(false);
template<> bool DefIVal<float>::isInit(false);
template<> bool DefIVal<double>::isInit(false);

template<> char DefIVal<char>::zeroVal(0);
template<> uchar DefIVal<uchar>::zeroVal(0);
template<> int DefIVal<int>::zeroVal(0);
template<> uint DefIVal<uint>::zeroVal(0);
template<> long DefIVal<long>::zeroVal(0);
template<> unsigned long DefIVal<unsigned long>::zeroVal(0);
template<> float DefIVal<float>::zeroVal(0.0);
template<> double DefIVal<double>::zeroVal(0.0);
