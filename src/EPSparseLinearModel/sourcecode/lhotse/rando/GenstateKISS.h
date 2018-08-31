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
 * Module: rando
 * Desc.:  Header class GenstateKISS
 * ------------------------------------------------------------------- */

#ifndef GENSTATEKISS_H
#define GENSTATEKISS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/GeneratorState.h"

//BEGINNS(rando)
  /**
   * Stores state for KISS generator 'GeneratorKISS'.
   * KISS uses three hidden chains of type unsigned int, so the state consists
   * of three unsigned int variables, called si,sj,sk here.
   * NOTE: For sk, the most significant bit is always 0.
   * <p>
   * NOTE: We use the assumptions (on the architecture) detailed in
   * 'FileUtils', esp. that unsigned int is 32 bits.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class GenstateKISS : public GeneratorState
  {
  protected:
    // Members

    uint si,sj,sk; // state variables

  public:
    // Public methods

    /**
     * Constructor.
     *
     * @param initSi Init. value for si
     * @param initSj Init. value for sj
     * @param initSk Init. value for sk. The most signif. bit is set to 0
     */
    GenstateKISS(uint initSi=0,uint initSj=0,uint initSk=0) : si(initSi),
    sj(initSj) {
      sk=initSk&0x7FFFFFFF;
    }

    /**
     * Copy constructor
     *
     * @param src Source object
     */
    GenstateKISS(const GenstateKISS& src) : si(src.si),sj(src.sj),sk(src.sk) {}

    /**
     * Draws a copy of state 'src' and stores it in this object.
     * ATTENTION: dynamic_cast broken for Matlab interface!
     *
     * @param src Source state
     */
    void copy(const GeneratorState* src) {
      const GenstateKISS* srcP=DYNCAST(const GenstateKISS,src);

      if (srcP==0)
	throw InvalidParameterException("'src' has wrong type");
      si=srcP->si; sj=srcP->sj; sk=srcP->sk;
    }

    void print() const {
      char msg[200];
      sprintf(msg,"PRN generator state: GenstateKISS\n  I=%u\n  J=%u\n  K=%u",
	      si,sj,sk);
      printMsgStdout(msg);
    }

    /**
     * The format is:
     *   <I-value>:<J-value>:<K-value>
     */
    void initFromString(const char* str) {
      uint tsi,tsj,tsk;
      if (sscanf(str,"%u:%u:%u",&tsi,&tsj,&tsk)!=3)
	throw InvalidParameterException(EXCEPT_MSG(""));
      si=tsi; sj=tsj; sk=tsk&0x7FFFFFFF;
    }

    void getString(string& str) const {
      char sbuff[80];
      sprintf(sbuff,"%u:%u:%u",si,sj,sk);
      str=sbuff;
    }

    /**
     * Initialises state to new values
     *
     * @param initSi New value for si
     * @param initSj New value for sj
     * @param initSk New value for sk. The most signif. bit is set to 0
     */
    virtual void set(uint initSi,uint initSj,uint initSk) {
      si=initSi; sj=initSj; sk=initSk&0x7FFFFFFF;
    }

    /**
     * @param siP si ret. here
     * @param sjP sj ret. here
     * @param skP sk ret. here
     */
    virtual void get(uint* siP,uint* sjP,uint* skP) {
      *siP=si; *sjP=sj; *skP=sk;
    }

    uint stateI() const {
      return si;
    }

    void stateI(uint vsi) {
      si=vsi;
    }
    
    uint stateJ() const {
      return sj;
    }

    void stateJ(uint vsj) {
      sj=vsj;
    }
    
    uint stateK() const {
      return sk;
    }

    void stateK(uint vsk) {
      sk=vsk&0x7FFFFFFF;
    }

    /**
     * Uses output of 'time' to initialise the state. The value for si is
     * computed directly from this output, sj and sk are initialised to
     * values obtained from si by running the congruental generator of KISS
     * 3 and 4 times.
     * ==> For serious simulations, this initialisation is NOT appropriate!
     *     Especially, the values ret. by 'time' are rather small.
     */
    void initFromTime() {
      si=(uint) time(0);
      uint temp=69069*si+23606797;
      temp=69069*temp+23606797; sj=69069*temp+23606797;
      sk=(69069*sj+23606797)&0x7FFFFFFF;
    }
  };
//ENDNS

#endif
