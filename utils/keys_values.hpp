#ifndef KEYS_VALUES_HPP
#define KEYS_VALUES_HPP

#include <stdexcept>
#include <map>
#include <complex>
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
#include <boost/tuple/tuple.hpp>


namespace {

  using namespace std;
  typedef complex<double> CD;
  using std::transform;
  using boost::any;
  using boost::any_cast;
  using boost::lexical_cast;
  using boost::algorithm::is_any_of;
  using boost::algorithm::split;
  using boost::algorithm::trim_copy;
  using boost::algorithm::trim;
  using boost::algorithm::token_compress_on;
  using boost::tuple;
  using boost::make_tuple;
  using boost::get;
  
  typedef const string CS;
}

/**
 * convert string to type T.
 * If some types are given such as T,U,...
 * return tuple<T,U,...> value.
 */
template<class T>
T ConvertData(const string& str);
template<class T, class U>
tuple<T, U> ConvertData(CS& str, CS& sep );
template<class T, class U, class V>
tuple<T,U,V> ConvertData(CS& str, CS& sep);
template<class T, class U, class V, class W>
tuple<T,U,V,W> ConvertData(CS& str, CS& sep);

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
  // return value with type T.
  // if casting T failed, throw an exception.
  // If unsuported T is called, compile failed.
  template<class T> T Get(const string& k, int i) const;
  template<class T> T Get(const string& k) const;
  //
  // return number of values for key.
  int Count(const string&) const;
  //
  // return existance of key. this function does not throw
  // any exceptions.
  bool ExistKey(const string&) const;
  //
  // same behavior with ExistKey but, if key does not
  // exist, throw exception.
  void CheckExistKey(const string&) const;
  void CheckIndex(const string&, int) const;

  // --------- Setter ----------
  //
  // add key k with value t.
  template<class T>
  void Add(const string& k, T t);
  //
  // add key k with value v. v is converted to some type.
  void AddAtomConverting(CS& k, CS& v);
  //
  // Convert data type whose key is k.
  template<class T> void ConvertValues(CS&);
  template<class T, class U> void ConvertValues(CS&);
  template<class T, class U, class V> void ConvertValues(CS&);
  template<class T, class U, class V, class W> 
  void ConvertValues(CS&);
  //
  // set value if the key is not exist.
  // this will is use for default values
  template<class T>
  bool SetIfNull(CS& key, T val);
  //
  // read line and add key and value as string data
  void ReadLine(CS& line);
  //
  // read file using ifstream
  void Read(ifstream& ifs);
};

#include "keys_values_priv.hpp"

#endif
