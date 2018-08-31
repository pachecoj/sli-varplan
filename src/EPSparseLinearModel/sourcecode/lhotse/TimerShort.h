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
 * Desc.:  Header class TimerShort
 * ------------------------------------------------------------------- */

#ifndef TIMERSHORT_H
#define TIMERSHORT_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h"

/**
 * Implements simple stopwatch for time measurements for fairly brief
 * time intervals. The function 'clock' from time.h (ctime) is used,
 * which measures processer time used by the program, but not by its
 * children.
 * <p>
 * ATTENTION: This is probably not portable at all. Claims to be ANSI C,
 * works for Linux.
 * ATTENTION: There is a wrap-around about every 72 minutes on a 32 bit
 * system, so the timer should be used for brief intervals only. There
 * is no simple way to detect whether one (or more) wrap-arounds have
 * taken place.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class TimerShort
{
protected:
  // Members

  time_t startT; // Timepoint when the counter was started
  time_t stopT;  // Timepoint when the counter was stopped
  bool isRun;    // Is the timer running?

public:
  // Constructors

  /**
   * Default constructor. Initially, the counter is not running and
   * 'startT', 'stopT' have the same value 0.
   */
  TimerShort() : startT(0),stopT(0),isRun(false) {}

  /**
   * Start the timer. If it is already running, it is reset.
   */
  void start() {
    startT=clock(); isRun=true;
  }

  /**
   * Stop the timer. If it is not running, this method does nothing. The
   * elapsed processor time (in secs) between start and stop point is
   * returned (see header comment about wrap-arounds!).
   *
   * @return Elapsed time
   */
  double stop() {
    if (isRun) {
      stopT=clock(); isRun=false;
    }

    return ((double) (stopT-startT))/CLOCKS_PER_SEC;
  }

  /**
   * Returns the elapsed processor time (in secs.) since the start point if
   * the timer is running. If the timer is not running, we return the time
   * between start and stop point. The timer is not stopped.
   *
   * @return See above
   */
  double getDiff() const {
    time_t actT=isRun?clock():stopT;

    return ((double) (actT-startT))/CLOCKS_PER_SEC;
  }

  /**
   * @return Is the counter running?
   */
  bool isRunning() const {
    return isRun;
  }
};

#endif
