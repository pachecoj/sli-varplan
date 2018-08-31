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
 * Desc.:  Definition classes ArgBlockType, ArgBlock
 * ------------------------------------------------------------------- */

#include "lhotse/ArgBlock.h"
#include "lhotse/CommandParser.h"
#include "lhotse/matrix/BaseVector.h"
#include "lhotse/matrix/StVector.h"

USING(matrix);

// Static members

const int ArgBlockType::presReqUncond;
const int ArgBlockType::presOptUncond;
const int ArgBlockType::presReqIfDep;
const int ArgBlockType::presReqIfNotDep;

// Public methods

ArgBlockType::ArgBlockType(const my_string& nameP,int typeP,void* defValP,
			   void* rngIvP,int pcheckP,int depParP)
  : name(nameP),type(typeP),pcheck(pcheckP),depPar(depParP),defVal(0),
    rngIv(0)
{
  if (typeP<0 || typeP>CommandParser::typeLast)
    throw InvalidParameterException("Invalid type");
  if (pcheck<0 || pcheck>presReqIfNotDep)
    throw InvalidParameterException("pcheck");
  defVal=ArgBlock::copyObj(defValP,typeP); // copy def. value obj.
  if (rngIvP!=0) {
    switch (typeP) {
    case CommandParser::typeInt:
      rngIv=new Interval<int>(*((Interval<int>*) rngIvP));
      break;
    case CommandParser::typeDouble:
      rngIv=new Interval<double>(*((Interval<double>*) rngIvP));
      break;
    case CommandParser::typeString:
      rngIv=0;
      break;
    case CommandParser::typeBool:
      rngIv=0;
      break;
    case CommandParser::typeLong:
      rngIv=0;
      break;
    case CommandParser::typeListInt:
      rngIv=new Interval<int>(*((Interval<int>*) rngIvP));
      break;
    case CommandParser::typeListDouble:
      rngIv=new Interval<double>(*((Interval<double>*) rngIvP));
      break;
    }
  }
}

ArgBlockType::~ArgBlockType()
{
  ArgBlock::freeObj(defVal,type);
  if (rngIv!=0) {
    switch (type) {
    case CommandParser::typeInt:
    case CommandParser::typeListInt:
      delete ((Interval<int>*) rngIv);
      break;
    case CommandParser::typeDouble:
    case CommandParser::typeListDouble:
      delete ((Interval<double>*) rngIv);
      break;
    case CommandParser::typeLong:
      delete ((Interval<long>*) rngIv);
      break;
    }
  }
}

bool ArgBlockType::checkIntVal(const void* obj) const
{
  if (rngIv==0) return true;
  if (type==CommandParser::typeString || type==CommandParser::typeBool)
    return true;
  int res=0;
  switch (type) {
  case CommandParser::typeInt:
    res=((Interval<int>*) rngIv)->check(*((int*) obj));
    break;
  case CommandParser::typeDouble:
    res=((Interval<double>*) rngIv)->check(*((double*) obj));
    break;
  case CommandParser::typeLong:
    res=((Interval<long>*) rngIv)->check(*((long*) obj));
    break;
  case CommandParser::typeListInt:
    if (!((BaseVector<int>*) obj)->checkBounds(*((Interval<int>*) rngIv)))
      res=1;
    break;
  case CommandParser::typeListDouble:
    if (!((StVector*) obj)->checkBounds(*((Interval<double>*) rngIv)))
      res=1;
    break;
  };

  return (res==0);
}

void ArgBlock::init(const CommandParser& pars,const my_string& prefNam,
		    const my_string& suffNam,const ArgBlock* block)
{
  int i,sz;

  prefName=prefNam; suffName=suffNam;
  // Has description table been setup?
  if (descr.isZero())
    throw InternalException("ArgBlock::init: Description table has not been setup!");
  // In a first sweep, we run over parameter of "pure" 'pcheck'
  sz=descr.size();
  value.changeRep(sz);
  given.changeRep(sz);
  for (i=0; i<sz; i++) given[i]=false;
  for (i=0; i<sz; i++)
    if (isPure(i)) readValue(pars,i,block);
  // Second sweep over the rest
  for (i=0; i<sz; i++)
    if (!isPure(i)) readValue(pars,i,block);
}

ArgBlock::~ArgBlock()
{
  int i,sz=descr.size();

  for (i=0; i<sz; i++)
    if (given[i]) {
      //cout << "Destroy " << i << endl; // DEBUG
      freeObj(value[i],descr[i]->type);
    }
}

int ArgBlock::getInt(const my_string& name) const
{
  int num=findPar(name,CommandParser::typeInt);
  return *((int*) value[num]);
}

double ArgBlock::getDouble(const my_string& name) const
{
  int num=findPar(name,CommandParser::typeDouble);
  return *((double*) value[num]);
}

const my_string& ArgBlock::getString(const my_string& name) const
{
  int num=findPar(name,CommandParser::typeString);
  return *((my_string*) value[num]);
}

bool ArgBlock::getBool(const my_string& name) const
{
  int num=findPar(name,CommandParser::typeBool);
  return *((bool*) value[num]);
}

long ArgBlock::getLong(const my_string& name) const
{
  int num=findPar(name,CommandParser::typeLong);
  return *((long*) value[num]);
}

const BaseVector<int>& ArgBlock::getListInt(const my_string& name) const
{
  int num=findPar(name,CommandParser::typeListInt);
  return *((BaseVector<int>*) value[num]);
}

const StVector& ArgBlock::getListDouble(const my_string& name) const
{
  int num=findPar(name,CommandParser::typeListDouble);
  return *((StVector*) value[num]);
}

void ArgBlock::set(const my_string& name,const void* valP)
{
  int num=findPar(name);
  if (!descr[num]->checkIntVal(valP)) {
    my_string msg("Value invalid for parameter '");
    msg+=getFullName(num); msg+="'";
    throw ArgBlockException(msg.c_str());
  }
  checkNewValue(valP,num);
  if (given[num])
    freeObj(value[num],descr[num]->type);
  value[num]=copyObj(valP,descr[num]->type);
  given[num]=true;
}

// Internal methods

int ArgBlock::findPar(const my_string& name,int typ) const
{
  LstIter it=nameToNum.find(name);
  if (it==nameToNum.end()) {
    my_string msg("Cannot find parameter '");
    msg+=name+"'";
    throw ArgBlockException(msg.c_str());
  }
  int num=it->second;
  if (typ!=-1 && typ!=descr[num]->type)
    throw ArgBlockException("Invalid return type");

  return num;
}

void ArgBlock::readValue(const CommandParser& pars,int num,
			 const ArgBlock* block)
{
  void* obj;
  const ArgBlockType* desc=descr[num];
  bool ok;
  my_string fname(getFullName(num));

  // Try to read
  try {
    obj=genObj(desc->type);
    pars.getValue(fname,desc->type,obj);
    given[num]=true;
  } catch (KeyNotFoundException ex) {
    given[num]=false;
  } catch (ParseException ex) {
    my_string msg("Value for parameter '");
    msg+=fname; msg+="' has invalid type";
    throw ArgBlockException(msg.c_str());
  }
  if (!given[num] && block!=0 && block->given[num]) {
    // Fill in from 'block'
    if (desc->type!=block->descr[num]->type)
      throw ArgBlockException("Cannot fill from 'block', has wrong type!");
    freeObj(obj,desc->type);
    obj=copyObj(block->value[num],desc->type); // draw copy
  }

  // Presence check
  if (!given[num] &&
      (desc->pcheck==ArgBlockType::presReqUncond ||
       (desc->pcheck==ArgBlockType::presReqIfDep && given[desc->depPar]) ||
       (desc->pcheck==ArgBlockType::presReqIfNotDep && !given[desc->depPar]))) {
    my_string msg("Parameter '");
    msg+=fname; msg+="' has to be given";
    throw ArgBlockException(msg.c_str());
  }
  // Interval check
  if (given[num] && !desc->checkIntVal(obj)) {
    my_string msg("Value for parameter '");
    msg+=fname; msg+="' violates range constraints";
    throw ArgBlockException(msg.c_str());
  }
  // OK, checks met
  if (given[num])
    value[num]=obj;
  else {
    freeObj(obj,desc->type);
    value[num]=desc->defVal; // def. value
  }
}

/*
 * If MATLAB_MEX is set, the Matlab MM has to be used even for the
 * elementary types, for which we cannot overload new/delete. They
 * have to be created in an awkward way.
 */
void* ArgBlock::genObj(int typ)
{
  switch (typ) {
#ifdef MATLAB_MEX
    void* newint,*newdouble,*newlong,*newbool;
#endif
  case CommandParser::typeInt:
#ifdef MATLAB_MEX
    newint=mxMalloc(sizeof(int));
    mexMakeMemoryPersistent(newint);
    int dumint; *((int*) newint)=dumint;
    return newint;
#else
    return new int();
#endif
  case CommandParser::typeDouble:
#ifdef MATLAB_MEX
    newdouble=mxMalloc(sizeof(double));
    mexMakeMemoryPersistent(newdouble);
    double dumdouble; *((double*) newdouble)=dumdouble;
    return newdouble;
#else
    return new double();
#endif
  case CommandParser::typeString:
    return new my_string();
  case CommandParser::typeLong:
#ifdef MATLAB_MEX
    newlong=mxMalloc(sizeof(long));
    mexMakeMemoryPersistent(newlong);
    int dumlong; *((long*) newlong)=dumlong;
    return newlong;
#else
    return new long();
#endif
  case CommandParser::typeBool:
#ifdef MATLAB_MEX
    newbool=mxMalloc(sizeof(bool));
    mexMakeMemoryPersistent(newbool);
    int dumbool; *((int*) newbool)=dumbool;
    return newbool;
#else
    return new bool();
#endif
  case CommandParser::typeListInt:
    return new BaseVector<int>();
  case CommandParser::typeListDouble:
    return new StVector();
  }
}

void* ArgBlock::copyObj(const void* obj,int typ)
{
  switch (typ) {
  case CommandParser::typeInt:
#ifdef MATLAB_MEX
    void* newint,*newdouble,*newlong,*newbool;
#endif
#ifdef MATLAB_MEX
    newint=mxMalloc(sizeof(int));
    mexMakeMemoryPersistent(newint);
    *((int*) newint)=*((int*) obj);
    return newint;
#else
    return new int(*((int*) obj));
#endif
  case CommandParser::typeDouble:
#ifdef MATLAB_MEX
    newdouble=mxMalloc(sizeof(double));
    mexMakeMemoryPersistent(newdouble);
    *((double*) newdouble)=*((double*) obj);
    return newdouble;
#else
    return new double(*((double*) obj));
#endif
  case CommandParser::typeString:
    return new my_string(*((my_string*) obj));
  case CommandParser::typeLong:
#ifdef MATLAB_MEX
    newlong=mxMalloc(sizeof(long));
    mexMakeMemoryPersistent(newlong);
    *((long*) newlong)=*((long*) obj);
    return newlong;
#else
    return new long(*((long*) obj));
#endif
  case CommandParser::typeBool:
#ifdef MATLAB_MEX
    newbool=mxMalloc(sizeof(bool));
    mexMakeMemoryPersistent(newbool);
    *((bool*) newbool)=*((bool*) obj);
    return newbool;
#else
    return new bool(*((bool*) obj));
#endif
  case CommandParser::typeListInt:
    return new BaseVector<int>(*((BaseVector<int>*) obj));
  case CommandParser::typeListDouble:
    return new StVector(*((StVector*) obj));
  }
}

void ArgBlock::freeObj(void* obj,int typ)
{
  if (obj!=0) {
    switch (typ) {
    case CommandParser::typeInt:
#ifdef MATLAB_MEX
      mxFree(obj);
#else
      delete ((int*) obj);
#endif
      break;
    case CommandParser::typeDouble:
#ifdef MATLAB_MEX
      mxFree(obj);
#else
      delete ((double*) obj);
#endif
      break;
    case CommandParser::typeString:
      delete ((my_string*) obj);
      break;
    case CommandParser::typeLong:
#ifdef MATLAB_MEX
      mxFree(obj);
#else
      delete ((long*) obj);
#endif
      break;
    case CommandParser::typeBool:
#ifdef MATLAB_MEX
      mxFree(obj);
#else
     delete ((bool*) obj);
#endif
      break;
    case CommandParser::typeListInt:
      delete ((BaseVector<int>*) obj);
      break;
    case CommandParser::typeListDouble:
      delete ((StVector*) obj);
      break;
    }
  } else
    cout << "WARNING: ArgBlock::freeObj: Called with zero pointer." << endl;
}

bool ArgBlock::presConsistent() const
{
  int i,sz=descr.size(),dep;

  for (i=0; i<sz; i++)
    if (!isPure(i)) {
      dep=descr[i]->depPar;
      if (dep<0 || dep>=sz || !isPure(dep)) return false;
    }
  return true;
}

void ArgBlock::checkGiven(int num) const
{
  if (!wasGiven(descr[num]->name)) {
    my_string msg("Parameter '");
    msg+=getFullName(num); msg+="' missing";
    throw ArgBlockException(msg.c_str());
  }
}
