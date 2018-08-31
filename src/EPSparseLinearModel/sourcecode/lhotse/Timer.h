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
 * Desc.:  Header class Timer
 * ------------------------------------------------------------------- */

#ifndef TIMER_H
#define TIMER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h"

/**
 * Implements simple stopwatch for time measurements for rather
 * long-running processes. We use the C library function 'time'
 * which returns the number of seconds passed since some point in
 * the past. Thus, the timer runs on a resolution of seconds.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class Timer
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
  Timer() : startT(0),stopT(0),isRun(false) {}

  /**
   * Start the timer. If it is already running, it is reset.
   */
  void start() {
    startT=time(0); isRun=true;
  }

  /**
   * Stop the timer. If it is not running, this method does nothing. The
   * elapsed time (in secs) between start and stop point is returned.
   *
   * @return Elapsed time
   */
  int stop() {
    if (isRun) {
      stopT=time(0); isRun=false;
    }

    return (int) (stopT-startT);
  }

  /**
   * Returns the elapsed time since the start point if the timer is running.
   * If the timer is not running, we return the time between start and stop
   * point. The timer is not stopped.
   *
   * @return See above
   */
  int getDiff() const {
    time_t actT=isRun?time(0):stopT;

    return (int) (actT-startT);
  }

  /**
   * @return Is the counter running?
   */
  bool isRunning() const {
    return isRun;
  }
};

#endif
