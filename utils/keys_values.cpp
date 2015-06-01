#include "keys_values.hpp"

bool IsInteger(const string& str) {
  
  bool acc = true;
  for(string::const_iterator it = str.begin();
      it != str.end(); ++it)
    acc = acc & isdigit(*it);
  
  return acc;
}
bool IsNumber(const string& str) {

  bool acc = true;
  for(string::const_iterator it = str.begin();
      it != str.end(); ++it)
    {
      char c = *it;
      acc = acc & (isdigit(c) ||c=='-'||c=='+'||c=='.'||c=='e');
    }
  return acc;
}

template<>
double KeysValues::Get<double>(CS& key, int i) const {
  CheckIndex(key, i);
  any x = dict_.find(key)->second[i];
  return any_cast<double>(x);
}
template<>
int KeysValues::Get<int>(CS& key, int i) const {

  this->CheckIndex(key, i);
  any x = dict_.find(key)->second[i];
  return any_cast<int>(x);
}
template<>
string KeysValues::Get<string>(CS& key, int i) const {

  /*
  ConstIt it = key_val_map_.find(key);
  if (it == key_val_map_.end())
    {
      string msg;
      msg += "Error: KeyVal::GetAsString\n";
      msg += "invalid key =" + key;
      throw msg;
    }

  any val = it->second;  
  std::type_info const & t = val.type();
  string res;

  if(t == typeid(int))
    res = lexical_cast<string>(any_cast<int>(val));
  else if(t == typeid(double))
    res = lexical_cast<string>(any_cast<double>(val));
  else if(t == typeid(string))
    res = any_cast<string>(val);
  else
    {
      string msg;
      msg += "Error: KeyVal::GetAsString\n";
      msg += "invalid type.";
      msg += "key is " + key;
    }
  return res;
  */
  return "this is not implemented";
}
int KeysValues::Count(CS& key) const {
  
  if(this->ExistKey(key))
    return dict_.find(key)->second.size();
  else
    return 0;
}
bool KeysValues::ExistKey (CS& key) const {
  return dict_.find(key) != dict_.end();
}
void KeysValues::CheckExistKey(CS& key) const {
  
  bool exist = ExistKey(key);
  if(not exist)
    {
      string msg ;
      msg += "Error: Key does not exi\n";
      msg += "key is " + key;
      throw invalid_argument(msg);
    }
}
void KeysValues::CheckIndex(CS& key, int i ) const {

  CheckExistKey(key);
  int num = this->Count(key);
  
  if(i > num - 1) {
    string msg;
    SUB_LOCATION(msg);
    msg += "Error: index i is out of range.";
    throw out_of_range(msg);
  }
  
}

void KeysValues::AddOld(CS& k, int t) {

  this->Add<int>(k, t);
  
  //  key_val_map_[k]=int(t); 
}
void KeysValues::AddOld(CS& k, double t) {
  
  this->Add<double>(k, t);
  //  key_val_map_[k]=double(t); 
}
void KeysValues::AddOld(CS& k, CS& t) {
  this->Add<string>(k, t);
  //  key_val_map_[k]=string(t); 
}

void KeysValues::AddCasting(CS& k, CS& v) {

  if(IsInteger(v))
    {
      this->Add(k, lexical_cast<int>(v));
    }
  else if(IsNumber(v))
    {
      this->Add(k, lexical_cast<double>(v));
    }
  else
    {
      this->Add(k, v);
    }
}
void KeysValues::AddSeparating(CS& k, CS& v, CS& sep) {

  if(v.find(sep) == string::npos) {
    
    // failed to find sep in v.
    // v is treated atomic value.
    this->AddCasting(k, v);
    
  } else {
    // find sep in v.
    // v is treated as tuple.
    
    this->Add();
  }
  
}

void KeyVal::Read(ifstream& ifs, string sep) {
  
  string line;
  while(getline(ifs, line))
    {
      if(line == "")
	continue;

      vector<string> str_vec;
      split(str_vec, line, is_any_of(sep));
      
      if(str_vec.size() != 2)
	{
	  string msg;
	  msg += "Error: KeyVal::Read\n";
	  msg += "error line is \n";
	  msg += line;
	  throw msg;
	}

      string key = trim_copy(str_vec[0]);
      string val = trim_copy(str_vec[1]);

      try
	{
	  this->AddCasting(key, val);
	}
      catch(...)
	{
	  string msg2;
	  msg2 += "Error: KeyVal::Read\n";
	  msg2 += "bad cast error\n";
	  msg2 += "error line is\n";
	  msg2 += line;
	  throw msg2;
	}
    }
}
/*
void KeyVal::Write(ofstream& ofs) const {
  
  for(ConstIt it = key_val_map_.begin();
      it != key_val_map_.end();
      ++it)
    {
      string k = it->first;
      string v = this->Get<string>(k);
      ofs << k << ": " << v << endl;
    }
}
*/
