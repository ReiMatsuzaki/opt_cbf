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

/*
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
CD KeysValues::Get<CD>(CS& k, int i) const {

  this->CheckIndex(k, i);
  any x = dict_.find(k)->second[i];
  return any_cast<CD>(x);
}
template<>
string KeysValues::Get<string>(CS& key, int i) const {

  this->CheckIndex(key, i);

  any val = dict_.find(key)->second[i];
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
}
*/
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
void KeysValues::AddAtomConverting(CS& k, CS& v) {

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
void KeysValues::ReadLine(CS& line) {

    // separate with sep_kv
    vector<string> kv;
    split(kv, line, is_any_of(sep_kv_));

    // check the structure
    if(kv.size() != 2) {
      string msg;
      SUB_LOCATION(msg);
      msg += "line must be 'key:value'";
      throw runtime_error(msg);
    }

    // copy
    string key = trim_copy(kv[0]);
    string val = trim_copy(kv[1]);
    // string key = kv[0];
    // string val = kv[1];

    // add key value as string
    this->Add<string>(key, val);
    
  }
void KeysValues::Read(ifstream& ifs) {
  
  string line;
  while(getline(ifs, line)) {

    if(line == "")
      continue;

    this->ReadLine(line);
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
