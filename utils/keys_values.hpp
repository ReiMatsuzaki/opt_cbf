#ifndef KEYS_VALUES_HPP
#define KEYS_VALUES_HPP

#include <stdexcept>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <algorithm>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include "macros.hpp"

namespace {

  using namespace std;
  using boost::any;
  using boost::any_cast;
  using boost::lexical_cast;
  using boost::tuple;
  using boost::make_tuple;
  typedef const string CS;
}

// class which inform bad template parameter
class BadTemplateParameter;

// data
class KeysValues {

private:
  // --------- Field -----------
  map<string, vector<any> > dict_;

public:
  // --------- Constructor ----
  KeysValues() {}
  // --------- Getter ---------
  //
  // Get:
  // return value with type T.
  // if casting T failed, throw an exception.
  // If unsuported T is called, compile failed.
  //
  // ExistKey:
  // return existance of key. this function does not throw
  // any exceptions.
  //
  // CheckExistKey:
  // same behavior with ExistKey but, if key does not
  // exist, throw exception.
  template<class T> T Get(const string& k, int) const {

    BadTemplateParameter a();
    /*
    string msg;
    SUB_LOCATION(msg);
    msg += "This type is not supported.";
    throw bad_typeid
    */
    return T();
  }
  int Count(const string&) const;
  bool ExistKey(const string&) const;
  void CheckExistKey(const string&) const;
  void CheckIndex(const string&, int) const;

  // --------- Setter ----------
  // Add:
  // 
  void AddOld(const string& k, int t);
  void AddOld(const string& k, double t);
  void AddOld(const string& k, const string& t);
  template<class T>
  void Add(const string& k, T t) {

    if(this->ExistKey(k)) {

      if(dict_[k][0].type() != typeid(T) ) {
	string msg;
	SUB_LOCATION(msg);
	msg += "\ntype mismatch. Inserting value type is ";
	msg += typeid(T).name();
	msg += "\n";
	throw invalid_argument(msg);	
      }

      dict_[k].push_back(t);
	
      
    } else {
      
      vector<any> ts;
      ts.push_back(t);
      dict_[k] = ts;
    }
  }
  void AddCasting(CS& k, CS& v);
  void AddSeparating(CS& k, CS& v, CS& sep);
};


#endif
