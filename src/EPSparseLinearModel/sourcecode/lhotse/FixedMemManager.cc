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
 * Desc.:  Definition class FixedMemManager
 * ------------------------------------------------------------------- */

#include "lhotse/FixedMemManager.h"

// Do sanity checks?
#define DO_SANITY_CHECKS

// Static members

const int FixedMemManager::maxSize;
const int FixedMemManager::lalMax;
const int FixedMemManager::fstPageSize;
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
const int FixedMemManager::slMax;
#endif

bool FixedMemManager::isInit(false);
bool FixedMemManager::doCleanup(false);
MAP_TYPE(int,PageLstElem) FixedMemManager::pageLst;
void** FixedMemManager::lastAllocLst(0);
int FixedMemManager::lalNum(0);
FixedMemManagerPage** FixedMemManager::laPageLst(0);
MAP_TYPE(int,bool) FixedMemManager::addSizeLst;
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
int* FixedMemManager::statLSize(0);
long* FixedMemManager::statLUsed(0);
int FixedMemManager::slNum(0);
long FixedMemManager::statNFast(0);
long FixedMemManager::statNSlow(0);
long FixedMemManager::statCpAlloc(0);
#endif

// Public static methods

void FixedMemManager::init()
{
  if (!isInit) {
    // Sanity check
    lalNum=0;
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
    slNum=0; statNFast=statNSlow=0;
    statCpAlloc=0;
#endif
#ifndef MATLAB_MEX
    lastAllocLst=(void**) malloc(lalMax*sizeof(void*));
    laPageLst=(FixedMemManagerPage**) malloc(lalMax*sizeof(FixedMemManagerPage*));
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
    statLSize=(int*) malloc(slMax*sizeof(int));
    statLUsed=(long*) malloc(slMax*sizeof(long));
#endif
#else
    void* ptr=mxMalloc(lalMax*sizeof(void*));
    mexMakeMemoryPersistent(ptr);
    lastAllocLst=(void**) ptr;
    ptr=mxMalloc(lalMax*sizeof(FixedMemManagerPage*));
    mexMakeMemoryPersistent(ptr);
    laPageLst=(FixedMemManagerPage**) ptr;
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
    ptr=mxMalloc(slMax*sizeof(int));
    mexMakeMemoryPersistent(ptr);
    statLSize=(int*) ptr;
    ptr=mxMalloc(slMax*sizeof(long));
    mexMakeMemoryPersistent(ptr);
    statLUsed=(long*) ptr;
#endif
#endif
    /* Sizes above 4*'maxSize' which should be represented by a page
       nevertheless */
    addSizeLst.insert(pair<int,bool>(104,true));

    isInit=true; // OK, initialized (this must be the last instruction!)
  }
}

void* FixedMemManager::alloc(int size)
{
  if (!isInit || doCleanup) {
    // MM not initialized (or cleanup phase): use standard MM
#ifndef MATLAB_MEX
    return malloc(size);
#else
    void* ptr=mxMalloc(size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
#endif
  }
  isInit=false; // switch off temporarily
  FixedMemSemaphore control; // controls flag

  int i;
  if (size<=0) throw InvalidParameterException(EXCEPT_MSG(""));
  // Make sure 'size'>=sizeof(FreeListElem), div. by 4
  if (size<sizeof(FreeListElem)) size=sizeof(FreeListElem);
  if ((i=size%4)!=0) size+=(4-i);
  bool doRep=(size<=maxSize*4);
  if (!doRep) {
    // Look in 'addSizeLst'
    MAP_CONSTITER(int,bool) it=addSizeLst.find(size);
    doRep=(it!=addSizeLst.end());
  }
  if (!doRep) {
    // Use standard MM
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
    // DEBUG: statistics
    for (i=0; i<slNum; i++) {
      if (statLSize[i]==size) {
	statLUsed[i]++;
	break;
      } else if (statLSize[i]>size) {
	memmove((void*) (statLSize+(i+1)),(void*) (statLSize+i),
		(slMax-1-i)*sizeof(int));
	memmove((void*) (statLUsed+(i+1)),(void*) (statLUsed+i),
		(slMax-1-i)*sizeof(long));
	statLSize[i]=size; statLUsed[i]=1;
	break;
      }
    }
    if (i==slNum && slNum<slMax) {
      statLSize[slNum]=size; statLUsed[slNum++]=1;
    }
#endif
#ifndef MATLAB_MEX
    return malloc(size);
#else
    void* ptr=mxMalloc(size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
#endif
  } else {
    MAP_ITER(int,PageLstElem) it=pageLst.find(size);
    bool newFst=(it==pageLst.end());
    FreeListElem* freeEl=0;
    if (!newFst) freeEl=(*it).second.fstFree;
    if (freeEl==0) {
      // Free list empty: need to allocate new page
      int pageNum=fstPageSize;
      if (!newFst) {
	// Use twice the size of last page
	pageNum=(*it).second.lstPage->num*2;
      }
      FixedMemManagerPage* newPage=
	allocPage(size,pageNum,newFst?0:((*it).second.lstPage),freeEl);
      if (newFst) {
	// Update 'pageLst': insert new entry for 'size'
	PageLstElem newEl;
	newEl.lstPage=newPage;
	newEl.fstFree=freeEl;
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
	newEl.statUsed=0;
#endif
	pair<MAP_ITER(int,PageLstElem),bool> ret=
	  pageLst.insert(pair<int,PageLstElem>(size,newEl));
	it=ret.first;
      } else {
	// Update 'pageLst': link in new free list
	(*it).second.lstPage=newPage;
	(*it).second.fstFree=freeEl;
      }
    }
    // Now, 'freeEl' points to a free element, and 'it' ref. to the entry
    // in 'pageLst' for 'size'
    (*it).second.fstFree=freeEl->next; // update free list
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
    (*it).second.statUsed++; // statistics (DEBUG)
#endif
    // Sanity check
    if (--(freeEl->page->numFree)<0)
      throw InternalException(EXCEPT_MSG(""));
    // Update 'lastAllocLst'
    if (lalNum==lalMax) {
      // Shift list
      memmove((void*) lastAllocLst,(void*) (lastAllocLst+1),
	      (lalMax-1)*sizeof(void*));
      memmove((void*) laPageLst,(void*) (laPageLst+1),
	      (lalMax-1)*sizeof(FixedMemManagerPage*));
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
      statCpAlloc++;
#endif
      lalNum--;
    }
    void* ptr=(void*) freeEl;
    lastAllocLst[lalNum]=ptr;
    laPageLst[lalNum++]=freeEl->page;
    return ptr;
  }
}

void FixedMemManager::dealloc(void* ptr)
{
  if (!isInit) {
    // MM not initialized: use standard MM
#ifndef MATLAB_MEX
    free(ptr);
#else
    mxFree(ptr);
#endif
    return;
  }
  isInit=false; // switch off temporarily
  FixedMemSemaphore control; // controls flag

  int i;
  FixedMemManagerPage* page=0; // page for 'ptr'
  // Search through 'lastAllocLst' (linear!)
  void** actPtr;
  for (i=lalNum-1,actPtr=lastAllocLst+(lalNum-1); i>=0 && *actPtr!=ptr;
       i--,actPtr--);
  if (i>=0) {
    // Found 'ptr' (fast dealloc): remove from lists (overwrite with last)
    page=laPageLst[i];
    *actPtr=lastAllocLst[--lalNum];
    laPageLst[i]=laPageLst[lalNum];
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
    statNFast++; // found it fast
#endif
  }

  // Search through all pages of all sizes
  MAP_ITER(int,PageLstElem) pageIt=pageLst.end(); // iter. into 'pageLst'
  if (page==0) {
    pageIt=pageLst.begin();
    while (page==0 && pageIt!=pageLst.end()) {
      FixedMemManagerPage* actPg=(*pageIt).second.lstPage;
      while (actPg!=0) {
	if (ptr>=actPg->buff && ptr<=actPg->lastPtr) {
	  page=actPg;
#ifdef DO_SANITY_CHECKS
	  // Sanity check
	  long li=((long) ptr)-((long) actPg->buff);
	  if (li%(actPg->size)!=0)
	    throw InvalidParameterException(EXCEPT_MSG("ptr"));
#endif
	  break;
	} else
	  actPg=actPg->pred;
      }
      ++pageIt;
    }
#ifdef FIXEDMEMMANAGER_COLLECT_STATS
    if (page!=0) statNSlow++; // found it slowly
#endif
  }

  // De-allocation.
  if (page==0) {
    // Not one of our entries: Use standard MM
#ifndef MATLAB_MEX
    free(ptr);
#else
    mxFree(ptr);
#endif
  } else {
    // 'page' contains the memory page
    // Sanity check
    if (++(page->numFree)>page->num)
      throw InternalException(EXCEPT_MSG(""));
    if (doCleanup && page->numFree==page->num)
      // Cleanup phase: dellocate page
      deallocPage(page,pageIt);
    else {
      // Lookup entry in 'pageLst' (if not already done)
      if (pageIt==pageLst.end()) {
	pageIt=pageLst.find(page->size);
#ifdef DO_SANITY_CHECKS
	// Sanity check
	if (pageIt==pageLst.end()) throw InternalException(EXCEPT_MSG(""));
#endif
      }
      PageLstElem* plEnt=&(*pageIt).second; // entry in 'pageLst'
      FreeListElem* freeEl=(FreeListElem*) ptr;
      freeEl->next=plEnt->fstFree; // link into free list (as first elem.)
      freeEl->page=page;
      plEnt->fstFree=freeEl;
    }
  }
}

#ifdef FIXEDMEMMANAGER_COLLECT_STATS
void FixedMemManager::printStats()
{
  char msg[100];
  long total,totalUsed;

  if (!isInit)
    throw WrongStatusException(EXCEPT_MSG("MM not initialized")); 
  isInit=false; // switch off temporarily
  FixedMemSemaphore control; // controls flag
  MAP_CONSTITER(int,PageLstElem) it=pageLst.begin();
  printMsgStdout("*** FixedMemManager usage statistics:");
  while (it!=pageLst.end()) {
    FixedMemManagerPage* page=(*it).second.lstPage;
    total=totalUsed=0;
    while (page!=0) {
      total+=(long) page->num;
      totalUsed+=(long) page->numFree;
      page=page->pred;
    }
    int sz=(*it).second.lstPage->size;
    totalUsed=(total-totalUsed)*((long) sz);
    total*=(long) sz;
    sprintf(msg,"  sz=%d: usage=%ld, total buff=%ld, used=%ld",sz,
	    (*it).second.statUsed,total,totalUsed);
    printMsgStdout(msg);
    ++it;
  }
  printMsgStdout("  Larger stuff:");
  for (int i=0; i<slNum; i++) {
    sprintf(msg,"  sz=%d: usage=%ld",statLSize[i],statLUsed[i]);
    printMsgStdout(msg);
  }
  sprintf(msg,"  Dealloc, found fast: %ld\n  Dealloc, found slow: %ld",
	  statNFast,statNSlow);
  printMsgStdout(msg);
  sprintf(msg,"  Shift lastAllocLst (alloc): %ld",statCpAlloc);
  printMsgStdout(msg);
}
#endif

// Internal methods

FixedMemManagerPage* FixedMemManager::allocPage(int sizeP,int numP,
						FixedMemManagerPage* predP,
						FreeListElem*& freeL)
{
#ifndef MATLAB_MEX
  FixedMemManagerPage* page=
    (FixedMemManagerPage*) malloc(sizeof(FixedMemManagerPage));
  char* buffP=(char*) malloc(sizeP*numP);
#else
  FixedMemManagerPage* page=
    (FixedMemManagerPage*) mxMalloc(sizeof(FixedMemManagerPage));
  mexMakeMemoryPersistent((void*) page);
  char* buffP=(char*) mxMalloc(sizeP*numP);
  mexMakeMemoryPersistent((void*) buffP);
#endif
  page->buff=(void*) buffP;
  page->num=numP; page->succ=0; page->pred=predP;
  page->numFree=numP; page->size=sizeP;
  page->lastPtr=(void*) (buffP+(sizeP*(numP-1)));
  // Build free list
  FreeListElem* act,*next=(FreeListElem*) buffP;
  freeL=next;
  for (int i=1; i<numP; i++) {
    act=next;
    act->page=page;
    next=(FreeListElem*) (buffP+(i*sizeP));
    act->next=next;
  }
  next->page=page;
  next->next=0;
  // Link in
  if (predP!=0) predP->succ=page;

  return page;
}

void FixedMemManager::deallocPage(FixedMemManagerPage* page,
				  MAP_ITER(int,PageLstElem) pgIt)
{
  int i;

  // Is page empty?
  if (page->numFree<page->num)
    throw WrongStatusException(EXCEPT_MSG("Page not empty"));
  // Find entry in 'pageLst'
  if (pgIt==pageLst.end()) {
    pgIt=pageLst.find(page->size);
#ifdef DO_SANITY_CHECKS
    // Sanity check
    if (pgIt==pageLst.end()) throw InternalException(EXCEPT_MSG(""));
#endif
  }
  // Clean up free list
  PageLstElem* plEnt=&(*pgIt).second;
#ifdef DO_SANITY_CHECKS
  // Sanity check
  if (plEnt->lstPage->size!=page->size)
    throw InternalException(EXCEPT_MSG(""));
#endif
  FreeListElem* freeEl=plEnt->fstFree,*freePred;
  while (freeEl!=0 && freeEl->page==page) freeEl=freeEl->next;
  freePred=plEnt->fstFree=freeEl;
  if (freeEl!=0) {
    freeEl=freeEl->next;
    for (; freeEl!=0; freeEl=freeEl->next) {
      if (freeEl->page==page)
	freePred->next=freeEl->next;
      else
	freePred=freeEl;
    }
  }
  // Modify 'pageLst'
  if (plEnt->lstPage==page) {
#ifdef DO_SANITY_CHECKS
    // Sanity check
    if (page->succ!=0) throw InternalException(EXCEPT_MSG(""));
#endif
    if ((plEnt->lstPage=page->pred)==0) {
      // No more pages for this block size
      pageLst.erase(pgIt); // remove entry
    }
  }
  // Reconnect pred. and succ.
  if (page->pred!=0) page->pred->succ=page->succ;
  if (page->succ!=0) page->succ->pred=page->pred;
  // Deallocate page buffer
#ifndef MATLAB_MEX
  free(page->buff);
  free((void*) page);
#else
  mxFree(page->buff);
  mxFree((void*) page);
#endif
}
