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
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
//#include <boost/tuple/tuple.hpp>
#include <tr1/tuple>
#include "macros.hpp"

namespace {

  using namespace std;
  using boost::any;
  using boost::any_cast;
  using boost::lexical_cast;
  using boost::algorithm::is_any_of;
  using boost::algorithm::split;
  using boost::algorithm::trim_copy;
  using boost::algorithm::token_compress_on;
  //  using boost::tuple;
  //  using boost::make_tuple;
  using std::tr1::tuple;
  using std::tr1::make_tuple;
  using std::tr1::get;
  typedef std::complex<double> CD;
  typedef const string CS;
}

// remove specific element
template<class T>
void RemoveElement( const T& ele, vector<T>* vec) {

  typename vector<T>::iterator it = vec->begin();
  while(it != vec->end()) {
    if( *it == ele) 
      it = vec->erase(it);
    else
      ++it;
  }
}

/**
 * convert string to type T.
 * If some types are given such as T,U,...
 * return tuple<T,U,...> value.
 */
template<class T>
T ConvertData(const string& str) {
  return lexical_cast<T>(str);
}
template<class T, class U>
tuple<T, U> ConvertData(CS& str, CS& sep ) {

  // split string
  vector<string> str_vec;
  split(str_vec, str, is_any_of(sep), token_compress_on);

  // remove black
  //  RemoveElement<string>(" ", &str_vec);

  // check size
  if(str_vec.size() < 2) {
    string msg;
    SUB_LOCATION(msg);
    msg += "\nthe number of splitted strings are less than 2\n";
    throw invalid_argument(msg);
  }

  // trim data
  string vt = trim_copy(str_vec[0]);
  string vu = trim_copy(str_vec[1]);
  
  // convert each value
  tuple<T, U> res(ConvertData<T>(vt),
		  ConvertData<U>(vu));
  return res;
  
}
template<class T, class U, class V>
tuple<T,U,V> ConvertData(CS& str, CS& sep) {

  // split string
  vector<string> str_vec;
  split(str_vec, str, is_any_of(sep), token_compress_on);

  // remove black element
  //  RemoveElement<string>(" ", &str_vec);

  // check size
  if(str_vec.size() < 3) {
    string msg;
    SUB_LOCATION(msg);
    msg += "\nthe number of splitted strings are less than 3\n";
    throw invalid_argument(msg);
  }

  // trim
  string tt = trim_copy(str_vec[0]);
  string ut = trim_copy(str_vec[1]);
  string vt = trim_copy(str_vec[2]);

  /*
  // convert each value
  string str_t, str_u, str_v;
  str_t = ConvertData<T>(tt);
  str_u = ConvertData<U>(ut);
  str_v = ConvertData<V>(vt);
  */

  tuple<T, U, V> res(ConvertData<T>(tt),
		     ConvertData<U>(ut),
		     ConvertData<V>(vt));
  return res;
}

// class which inform bad template parameter
//template<class T> struct BadTemplateParameter;

// data
class KeysValues {

private:
  // --------- Field -----------
  map<string, vector<any> > dict_;
  string sep_kv_;  // separator for key value
  string sep_val_; // separator for values

public:
  // --------- Constructor ----
  KeysValues(CS& _sep_kv, CS& _sep_val) {

    sep_kv_  = _sep_kv;
    sep_val_ = _sep_val;
    
  }
  
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
  template<class T> T Get(const string& k, int i) const {

    this->CheckIndex(k, i);
    any x = dict_.find(k)->second[i];
    return any_cast<T>(x);
    
  }
  int Count(const string&) const;
  bool ExistKey(const string&) const;
  void CheckExistKey(const string&) const;
  void CheckIndex(const string&, int) const;

  // --------- Setter ----------
  //
  // add key k with value t.
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
  //
  // add key k with value v. v is converted to some type.
  void AddAtomConverting(CS& k, CS& v);
  //
  // Convert data type whose key is k.
  template<class T>
  void ConvertValues(CS& k) {

    this->CheckExistKey(k);

    // check original is string type.
    // notice that every element of dict_[k] is same type.
    if(dict_[k][0].type() != typeid(string)) {
      string msg;
      SUB_LOCATION(msg);
      msg += "\nConverting values must be string\n";
      throw std::runtime_error(msg);
    }

    // loop for values for key k;
    for(vector<any>::iterator it = dict_[k].begin();
	it != dict_[k].end(); ++it) {
      
      string str = any_cast<string>(*it);
      
      T converted;

      try {
	converted = ConvertData<T>(str);
      } catch (exception& e) {
	string msg = e.what();
	msg += "\n";
	SUB_LOCATION(msg);
	msg += "\nfailed to conversion.\n";
	msg += "key: " + k + "\n";
	throw runtime_error(msg);
      }
      *it = converted;
    }
    
  }
  template<class T, class U>
  void ConvertValues(CS& k) {

    this->CheckExistKey(k);

    // check original is string type.
    // notice that every element of dict_[k] is same type.
    if(dict_[k][0].type() != typeid(string)) {
      string msg;
      SUB_LOCATION(msg);
      msg += "\nConverting values must be string\n";
      throw std::runtime_error(msg);
    }

    // loop for values for key k;
    for(vector<any>::iterator it = dict_[k].begin();
	it != dict_[k].end(); ++it) {
      string str = any_cast<string>(*it);

      *it = ConvertData<T,U>(str, sep_val_);
    }
  }
  template<class T, class U, class V>
  void ConvertValues(CS& k) {

    this->CheckExistKey(k);

    // check original is string type.
    // notice that every element of dict_[k] is same type.
    if(dict_[k][0].type() != typeid(string)) {
      string msg;
      SUB_LOCATION(msg);
      msg += "\nConverting values must be string\n";
      throw std::runtime_error(msg);
    }

    // loop for values for key k;
    for(vector<any>::iterator it = dict_[k].begin();
	it != dict_[k].end(); ++it) {
      string str = any_cast<string>(*it);
      *it = ConvertData<T,U,V>(str, sep_val_);
    }
  }
  //
  // read line and add key and value as string data
  void ReadLine(CS& line);
  //
  // read file using ifstream
  void Read(ifstream& ifs);
};


#endif
