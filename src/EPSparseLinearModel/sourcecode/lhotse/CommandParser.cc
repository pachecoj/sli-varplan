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
 * Desc.:  Definition class CommandParser
 * ------------------------------------------------------------------- */

#include "lhotse/CommandParser.h"
#include "lhotse/MachDep.h"

#include "lhotse/matrix/BaseVector.h"

USING(matrix);

const int CommandParser::maxLineLength;
const int CommandParser::maxListEntries;
const int CommandParser::typeInt;
const int CommandParser::typeDouble;
const int CommandParser::typeString;
const int CommandParser::typeBool;
const int CommandParser::typeLong;
const int CommandParser::typeListInt;
const int CommandParser::typeListDouble;
const int CommandParser::typeLast;

CommandParser::CommandParser(const my_string& fname)
{
  ArrayHandle<char> line(maxLineLength+1);
  char* act,*end,*find,*temp;
  int lineNo=0,numEntries=0;
  my_string key,val;
  ifstream is;
  char text[45];

  FileUtils::openFileRead(fname.c_str(),is);
  printMsgStdout("Parsing control file...");
  while (is) {
    is.getline(line,maxLineLength); lineNo++;
    // Strip leading,trailing whitespace
    act=line;
    while (isspace(*act)) act++;
    end=line+strlen(line)-1;
    while (end>=act && isspace(*end)) end--;
    if (end<act) continue; // empty line
    *(++end)=0;
    // Check if comment
    if (*act=='#') continue; // OK, comment
    // Look for seperating '='
    if ((find=strchr(act,'='))==0) {
      sprintf(text,"Parse error in line %5d of command file.",lineNo);
      throw ParseException(text);
    }
    temp=find-1;
    while (temp>=act && isspace(*temp)) temp--;
    // Note that keywords must not be empty, but values can be.
    if (temp<act) {
      sprintf(text,"Parse error in line %5d of command file.",lineNo);
      throw ParseException(text);
    }
    *(temp+1)=0; // overwrites either WS or the '='
    key=act; // located keyword
    act=find+1;
    while (act<end && isspace(*act)) act++;
    val=act;
    keyValLst.insert(pair<my_string,my_string>(key,val));
    numEntries++;
  }
  CLOSESTR(is);
  sprintf(text,"Found %d entries.",numEntries);
  printMsgStdout(text);
}

void CommandParser::getValue(const my_string& key,int type,void* value,
			     const IntVal* ival) const
{
  int i,j,num,val;
  char* act,*begin;
  ArrayHandle<char*> startPos(maxListEntries);
  ArrayHandle<char> sbuff(maxLineLength+1);

  if (type<0 || type>typeLast)
    throw KeyNotFoundException("Invalid type code");
  LstIter it=keyValLst.find(key);
  if (it==keyValLst.end()) {
    my_string msg="Cannot find parameter '";
    msg+=key; msg+="'";
    throw KeyNotFoundException(msg.c_str());
  } else {
    const my_string& valStr=it->second;
    try {
      if (type==typeInt) {
	val=MachDep::stringToInt(valStr.c_str());
	const Interval<int>* iv=DYNCAST(const Interval<int>,ival);
	if (iv!=0 && iv->check(val)!=0) {
	  my_string msg="Parameter '";
	  msg+=key; msg+="' out of range";
	  throw OutOfRangeException(msg.c_str());
	}
	*((int*) value)=val;
      } else if (type==typeDouble) {
	double dval=MachDep::stringToDouble(valStr.c_str());
	const Interval<double>* iv=DYNCAST(const Interval<double>,ival);
	if (iv!=0 && iv->check(dval)!=0) {
	  my_string msg="Parameter '";
	  msg+=key; msg+="' out of range";
	  throw OutOfRangeException(msg.c_str());
	}
	*((double*) value)=dval;
      } else if (type==typeString) {
	*((my_string*) value)=valStr;
      } else if (type==typeBool) {
	val=MachDep::stringToInt(valStr.c_str());
	*((bool*) value)=(val!=0);
      } else if (type==typeLong) {
	long lval=MachDep::stringToLong(valStr.c_str());
	const Interval<long>* iv=DYNCAST(const Interval<long>,ival);
	if (iv!=0 && iv->check(lval)!=0) {
	  my_string msg="Parameter '";
	  msg+=key; msg+="' out of range";
	  throw OutOfRangeException(msg.c_str());
	}
	*((long*) value)=lval;
      } else if (type==typeListInt || type==typeListDouble) {
	// 'values[i]' contains the elements, sep. by whitespace. We first
	// parse the string to locate the starting points. Note that
	// 'values[i]' cannot have leading or trailing WS.
	strcpy(sbuff.p(),valStr.c_str());
	num=0; act=begin=sbuff;
	while (*act!=0) {
	  while (!isspace(*act) && *act!=0) act++;
	  startPos[num++]=begin;
	  if (*act!=0) {
	    *act=0; // mark end
	    act++;
	    while (isspace(*act)) act++;
	    begin=act;
	  }
	}
	if (type==typeListInt) {
	  BaseVector<int>* vec=(BaseVector<int>*) value;
	  vec->fill(num,0);
	  const Interval<int>* iv=DYNCAST(const Interval<int>,ival);
	  for (j=0; j<num; j++) {
	    (*vec)[j]=val=MachDep::stringToInt(startPos[j]);
	    if (iv!=0 && iv->check(val)!=0) {
	      my_string msg="Parameter '";
	      msg+=key; msg+="' out of range";
	      throw OutOfRangeException(msg.c_str());
	    }
	  }
	} else {
	  BaseVector<double>* vec=(BaseVector<double>*) value;
	  vec->fill(num,0.0);
	  const Interval<double>* iv=DYNCAST(const Interval<double>,ival);
	  double dval;
	  for (j=0; j<num; j++) {
	    (*vec)[j]=dval=MachDep::stringToDouble(startPos[j]);
	    if (iv!=0 && iv->check(dval)!=0) {
	      my_string msg="Parameter '";
	      msg+=key; msg+="' out of range";
	      throw OutOfRangeException(msg.c_str());
	    }
	  }
	}
      }
    } catch (ParseException ex) {
      my_string msg="Parse exception for value of parameter '";
      msg+=key; msg+="'";
      throw ParseException(msg.c_str());
    }
  }
}
