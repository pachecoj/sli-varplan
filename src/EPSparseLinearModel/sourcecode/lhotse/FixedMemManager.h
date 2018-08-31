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
 * Desc.:  Header class FixedMemManager
 * ------------------------------------------------------------------- */

#ifndef FIXEDMEMMANAGER_H
#define FIXEDMEMMANAGER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/*
 * TODO:
 * - A global replacement of 'new', 'delete' seems to be critical, leads
 *   to segmentation fault at least when the program terminates.
 *   ==> try rather to identify CORE CLIENTS of this manager, such as
 *       'Handle', 'ArrayHandle', mask vectors/matrices (maybe all
 *       vectors/matrices). Then, replace new/delete for these only!
 * - Method for compressing the pages. In the moment, empty pages are
 *   not deallocated
 */

#include "lhotse/global.h"

// Collect statistics? This is a small overhead

//#define FIXEDMEMMANAGER_COLLECT_STATS

struct FixedMemManagerPage;
struct FreeListElem;
struct PageLstElem;

/**
 * Element of free list.
 * ATTENTION: The size of this type is the minimum block size.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
struct FreeListElem
{
  FreeListElem* next;        // next element, or 0
  FixedMemManagerPage* page; // page for elem.
};

/**
 * Internal struct for 'FixedMemManager'.
 * Manages a page for the MM. Pages are doubly linked and manage an entry
 * buffer, elements of size 'size'.
 * 'lastPtr' is pointing to the last element of the array, needed for
 * fast checks in 'FixedMemManager::dealloc'.
 * <p>
 * NOTE: Why is this not a class with constructor? The problem is that in
 * this case, the overloaded 'new' is used to alloc. memory for objects of
 * 'FixedMemManagerPage' which may lead to trouble!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
struct FixedMemManagerPage
{
  void* buff;                 // array
  int size;                   // element size
  int num;                    // number of elements in the array
  int numFree;                // number of free elements in array
  FixedMemManagerPage* succ;  // successor array; or 0
  FixedMemManagerPage* pred;  // predecessor array; or 0
  void* lastPtr;              // s.a.
};

/**
 * Value type for 'FixedMemManager::pageLst'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
struct PageLstElem
{
  FixedMemManagerPage* lstPage; // last page for block size
  FreeListElem* fstFree;        // first elem. of free list (or 0)
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
  long statUsed;                // how many alloc? (DEBUG)
#endif
};

/**
 * Efficient memory manager for blocks of identical size. The use of
 * this manager is more efficient than the normal heap, and there
 * is no overhead for organisation.
 * The block (or element) sizes are k*4 bytes, k <= 'maxSize'. The block
 * size is >= the size of 'FreeListElem'. See doc/fixedMemManager.txt for
 * details.
 * <p>
 * Initialisation:
 * Note that 'new', 'delete' are globally overwritten to call 'alloc' and
 * 'dealloc' here. This creates the problem that memory may be allocated
 * BEFORE the manager is ready, esp. within the initialization method
 * 'init'.
 * Therefore, 'alloc' and 'dealloc' check the 'isInit' flag and use the
 * standard MM if this is false. They do not call 'init'.
 * 'init' has to be called explicitly, and only from this point on this MM
 * will be used.
 * NOTE: 'dealloc' uses the standard MM for adresses it cannot account for,
 * so memory allocated before 'init' was called is dealloc. properly.
 * <p>
 * Semaphore:
 * The 'isInit' flag acts as a semaphore. It is reset at the beginning of
 * 'alloc' / 'dealloc', and set again at the end. Makes sure this MM is not
 * used for anything within these methods.
 * <p>
 * How it works:
 * 'pageLst' maps existing block sizes to 'PageLstElem' descriptors, cont.
 *  ref. to last page and to 1st elem. of free list (0 if no free elem.).
 * The free list elem. are stored in free page slots, using the type
 * 'FreeListElem' (cont. pointer to next free list elem., and pointer to
 * page for this elem.). All free elem. in all ex. pages are in the corr.
 * free list. If a free list is empty, 'alloc' creates a new page.
 * Pages for a blocksize are collected in a doubly linked list, nodes of
 * type 'FixedMemManagerPage'. These also store the number of free slots
 * on the page (which are listed in the free list). Pages are only ref. to
 * by 'pageLst' entries.
 * <p>
 * List of last recently allocated:
 * We keep this list to avoid expensive search through all pages of
 * all sizes for 'ptr' in 'dealloc'.
 * 'lastAllocLst' of size 'lalNum' (max. 'lalMax') contains ptr. for
 * the last recently allocated elements, 'laPageLst' provides the
 * corr. pages for these elements.
 * NOTE: These are NOT kept in sorted order and 'dealloc' searches
 * through 'lastAllocLst' linearly: only OK if 'lalMax' not large!
 * Once 'lalNum'=='lalMax' in 'alloc', the arrays are shifted to the
 * left, the leftmost elements are removed.
 * If an element is found in 'dealloc', the current rightmost element
 * takes its place. Avoids copying there, but also jumbles the order:
 * elements alloc. more recently can be left of earlier ones.
 * ==> TODO: Timing stats on whether keeping a sorted list would
 *     be more efficient. If 'lalMax' is not too small, shifts in
 *     'alloc' do not occur very often.
 * <p>
 * Caveats:
 * This MM is to replace the standard one globally, by overloading
 * new/delete. The problem is that we have to make sure it is
 * destroyed AFTER all its users. Our solution is to make all members,
 * methods static. There is no destructor, so there is no clean-up
 * when the program finishes.
 * ==> Creating a global 'FixedMemManager' variable did not work,
 *     because it was destroyed before others, resulting in seg.
 *     faults at program termination. Also, if MATLAB_MEX is set,
 *     the MM object must be static.
 * <p>
 * ATTENTION: This is not a secure MM (focus is on simplicity).
 * For example, the following actions are not detected and lead to
 * unpredictable errors downstream:
 * - deallocating the same region more than once
 * - calling 'dealloc' with a pointer which falls in the buffer
 *   region of one of the existing pages
 * <p>
 * Element size distribution from some experiments:
 * Typical ivmevmax:
 * In a typical LHOTSE application, the most used element size is
 * 80 bytes (70%), which is probably a mask vector, followed by
 * 8 bytes (20%), probably a handle(?) Also a lot at 16, 32 bytes.
 * Therefore, 'maxSize'==20 seems useful. The next follow-up is
 * 104 bytes (1%), repres. in 'addSizeLst'. Then comes 200 bytes
 * (0.2%), not repres.
 * For all repres. sizes, only a single page is allocated (may be
 * different in different applications).
 * Typical IVMMULTICLASS with C++ RUNEPPHASE:
 * - 80  (58%)
 * - 8   (34%)
 * - 16  (3.3%)
 * - 32  (1.9%)
 * - 40  (1.2%)
 * - 104 (0.9%)
 * All other sizes below 0.02%. Only 0.003% of all dealloc. are slow.
 * <p>
 * Slow deallocations:
 * Have to keep "slow dealloc" (one has to search through all pages)
 * to a minimum. Strategies:
 * - make sure ALL frequently requested sizes are served by a page.
 *   All sizes up to 4*'maxSize' are served, furthermore the ones
 *   entered in 'addSizeLst'
 *   ==> Run your applic. with FIXEDMEMMANAGER_COLLECT_STATS def.,
 *       check whether frequent sizes are missed!
 *   NOTE: A 'dealloc' for a size which is not served, is always slow!
 *   NOTE: The statistics only track the 'slMax' smallest unrepres.
 *   sizes, a run with fairly large 'slMax' may be req. to find all
 *   frequently used sizes.
 * - increase 'lalMax' ==> Leads to more time req. to find entries in
 *   'lastAllocLst' (linear search in the moment)!
 *
 * The lastAllocLst list:
 * This is kept as linear vector in the moment, kind of FIFO. But the
 * removal in 'dealloc' replaces the element to be removed by the
 * rightmost one, instead of shifting all elements to the right. This
 * leads to a small increase in the number of "slow dealloc", because
 * recently alloc. elements move to the left and are kicked out too
 * fast.
 * ==> Do more timing tests! Try to replace 'lastAllocLst' by a sorted
 *     data structure. How to do FIFO then?? Compare using timing tests.
 * <p>
 * Deallocation:
 * During normal use, pages are not deallocated (TODO: method for
 * compression!). Once method 'cleanup' is called, the MM switches into
 * cleanup phase by setting the 'doCleanup' flag. In this phase (which
 * is not left anymore), calls to 'alloc' use the standard MM, and in
 * 'dealloc' empty pages are deallocated.
 * NOTE: There is no mechanism for cleaning up leftover stuff for which
 * 'dealloc' is not called.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class FixedMemManager
{
  friend class FixedMemSemaphore;

public:
  // Constants

  static const int maxSize    =20;
  static const int lalMax     =64;
  static const int fstPageSize=16;

protected:
  // Variables

  static bool isInit;                       // Members properly init.?
  static bool doCleanup;                    // Final cleanup phase?
  static MAP_TYPE(int,PageLstElem) pageLst; // list for ex. block sizes
  static void** lastAllocLst;               // s.a.
  static int lalNum;                        // size of 'lastAllocLst'
  static FixedMemManagerPage** laPageLst;   // s.a.
  static MAP_TYPE(int,bool) addSizeLst;     // s.a.

#ifdef FIXEDMEMMANAGER_COLLECT_STATS
  // DEBUG: Stats for sizes larger than 4*'maxSize':
  static int* statLSize;
  static long* statLUsed;
  static int slNum;
  static const int slMax=16;
  static long statNFast;     // Dealloc, found in 'lastAllocLst'
  static long statNSlow;     // Dealloc, had to look through all pages
  static long statCpAlloc;   // Shift 'lastAllocLst' in 'alloc'
#endif

public:
  /**
   * Has to be called in order to kick off the manager. Before 'init' is
   * called, the standard MM is used.
   * If 'isInit'==false, all static members are initialized properly.
   * 'isInit' is set to true.
   */
  static void init();

  /**
   * To be called once this MM is not used anymore
  /**
   * @return Is MM initialized?
   */
  static bool isMMInit() {
    return isInit;
  }

  /**
   * Allocates element of size 'size' (bytes). In fact, the region returned
   * may be slightly larger.
   *
   * @param size Size in bytes
   * @return     Pointer to region
   */
  static void* alloc(int size);

  /**
   * Deallocates region at 'ptr'. We first search through 'lastAllocLst',
   * then through all pages of all sizes. If nothing is found, the standard
   * MM is used for dealloc.
   *
   * @param ptr Pointer to region
   */
  static void dealloc(void* ptr);

#ifdef FIXEDMEMMANAGER_COLLECT_STATS
  /**
   * Prints statistics on number of 'alloc' calls
   */
  static void printStats();
#endif

protected:
  // Internal methods

  /**
   * Allocates new page of type 'FixedMemManagerPage' and links it behind
   * 'predP' (succ. is 0). This involves linking all entries linearly in a
   * new free list.
   *
   * @param sizeP Element size
   * @param numP  Number of entries
   * @param predP Predecessor
   * @param freeL Pointer to first elem. of new free list ret. here
   * @return      Pointer to new page
   */
  static FixedMemManagerPage* allocPage(int sizeP,int numP,
					FixedMemManagerPage* predP,
					FreeListElem*& freeL);

  /**
   * Deallocates memory page. The cleanup comprises:
   * - check 'numFree', must be the same as 'num'. Exception otherwise
   * - run through free list for the block size, remove all entries for
   *   the page. Needs linear time!
   * - reconnect pred. and succ. page
   * - modify 'pageLst' if this page is the last one
   *
   * @param page Page to dealloc
   * @param pgIt Optional. Iterator pointing to 'pageLst' entry for block
   *             size of page. Def.: pageLst.end()
   */
  static void deallocPage(FixedMemManagerPage* page,
			  MAP_ITER(int,PageLstElem) pgIt=pageLst.end());
};

/**
 * Helper class used within 'FixedMemManager' methods. Objects control the
 * semaphore 'FixedMemManager::isInit'. The destructor sets the flag.
 * By deleting the flag and constructing an object at the beginning of a
 * 'FixedMemManager' method, the flag is controlled automatically: it is
 * set again once the method is left.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class FixedMemSemaphore
{
  friend class FixedMemManager;

  FixedMemSemaphore() {}

  ~FixedMemSemaphore() {
    FixedMemManager::isInit=true;
  }
};

#endif
