#ifndef KEYS_VALUES_PRIV_HPP
#define KEYS_VALUES_PRIV_HPP

#include "macros.hpp"

template<class T>
T ConvertData(const string& str) {

  T res;
  try {
    res = lexical_cast<T>(str);
  } catch (exception& e) {
    string msg;
    SUB_LOCATION(msg);
    msg += "ConvertData failed.";
    msg += "Original string : ";
    msg += str;
    throw runtime_error(msg);
  }
  return res;
}

template<class T, class U>
tuple<T, U> ConvertData(CS& str, CS& sep ) {

  // split string
  vector<string> str_vec;
  split(str_vec, str, is_any_of(sep), token_compress_on);

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

template<class T> 
T KeysValues::Get(const string& k, int i) const {

  this->CheckIndex(k, i);
  any x = dict_.find(k)->second[i];
  T   t;

  try {
    t = any_cast<T>(x);
  } catch(exception& e) {
    string msg;
    SUB_LOCATION(msg);
    msg += "\nfailed to any_cast\n";
    msg += "key: ";
    msg += k;
    throw runtime_error(msg);
  }
  return t;
    
}
template<class T> 
T KeysValues::Get(const string& k) const {
  return this->Get<T>(k, 0);
}

template<class T>
void KeysValues::Add(const string& k, T t) {

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
template<class T>
bool KeysValues::SetIfNull(CS& key, T val) {
  
  if(this->ExistKey(key))
    return false;

  this->Add(key, val);

  return true;
}

template<class T>
void KeysValues::ConvertValues(CS& k) {

    this->CheckExistKey(k);

    // check original is string type.
    // notice that every element of dict_[k] is same type.
    bool is_string = dict_[k][0].type() == typeid(string);
    bool is_type   = dict_[k][0].type() == typeid(T);

    if(is_type)
      return;
    
    if(!is_string) {
      string msg;
      SUB_LOCATION(msg);
      msg += "\nConverting values must be string \n";
      throw std::runtime_error(msg);
    }

    // loop for values for key k;
    for(vector<any>::iterator it = dict_[k].begin();
	it != dict_[k].end(); ++it) {

      string str;
      try {
	str = any_cast<string>(*it);
      } catch(exception& e) {
	string msg;
	SUB_LOCATION(msg);
	msg += "\nfailed to cast\n";
	msg += "key: ";
	msg += k;
	throw runtime_error(msg);
      }
      
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
void KeysValues::ConvertValues(CS& k) {

    this->CheckExistKey(k);
    // check original is string type.
    // notice that every element of dict_[k] is same type.
    bool is_string = dict_[k][0].type() == typeid(string);
    bool is_type   = dict_[k][0].type() == typeid(tuple<T,U>);

    if(is_type) {
      return;
    }

    if(!is_string) {
      string msg;
      SUB_LOCATION(msg);
      msg += "\nConverting values must be string \n";
      throw std::runtime_error(msg);
    }

    // loop for values for key k;
    for(vector<any>::iterator it = dict_[k].begin();
	it != dict_[k].end(); ++it) {
      string str = any_cast<string>(*it);

      try {
	*it = ConvertData<T,U>(str, sep_val_);
      } catch (exception& e) {
	string msg;
	SUB_LOCATION(msg);
	msg += "\n failed to ConvertData. ";
	msg += "key: ";
	msg += k;
	throw runtime_error(msg);
      }
    }
  }
template<class T, class U, class V>
void KeysValues::ConvertValues(CS& k) {

    this->CheckExistKey(k);

    // check original is string type.
    // notice that every element of dict_[k] is same type.
    bool is_string = dict_[k][0].type() == typeid(string);
    bool is_type   = dict_[k][0].type() == typeid(tuple<T,U,V>);

    if(is_type) {
      return;
    }

    if(!is_string) {
      string msg;
      SUB_LOCATION(msg);
      msg += "\nConverting values must be string\n";
      throw std::runtime_error(msg);
    }

    // loop for values for key k;
    for(vector<any>::iterator it = dict_[k].begin();
	it != dict_[k].end(); ++it) {
      string str = any_cast<string>(*it);

      try {
	*it = ConvertData<T,U,V>(str, sep_val_);
      } catch (exception& e) {
	string msg;
	SUB_LOCATION(msg);
	msg += "\n failed to ConvertData. ";
	msg += "key: ";
	msg += k;
	throw runtime_error(msg);
      }      
    }
  }

#endif
