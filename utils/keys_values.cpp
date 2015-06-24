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

KeysValues::KeysValues(CS& _sep_kv, CS& _sep_val) {

  sep_kv_  = _sep_kv;
  sep_val_ = _sep_val;
    
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

void KeysValues::AddNull(CS& k)  {
    
  if(not this->ExistKey(k)) {
    vector<any> ts;
    dict_[k] = ts;
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

template<>
bool ConvertData<bool>(const string& str) {

  string tmp = str;
  trim(tmp);
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  
  if(tmp == "TRUE" || tmp == "T") 
    return true;
  else if(tmp == "FALSE" || tmp == "F")
    return false;

  string msg;
  SUB_LOCATION(msg);
  msg += "\nConvertData<bool> failed. ";
  msg += "\nOriginal string : ";
  msg += str;
  msg += "\nTransformed string : \"";
  msg += tmp;
  msg += "\"\n";
  throw runtime_error(msg);
  
}


