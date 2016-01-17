// Copyright 2012 Victor Bittorf, Chris Re
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include <string>
#include <iostream>
#include <iomanip>

#include "frontend_util.h"

namespace hazy {
namespace hogwild {

option* convert_extended_options(const extended_option* exs)  {
  using namespace std;
  int nOptions = 0; 
  while(exs[nOptions++].name != NULL);
  struct option* opts =  new option[nOptions];
  for(int i = 0; i < nOptions; i++) {
    opts[i].name    = exs[i].name;
    opts[i].has_arg = exs[i].has_arg;
    opts[i].flag    = exs[i].flag;
    opts[i].val     = exs[i].val;
  }
  return opts;
}

void print_usage(const extended_option *exs, char* sysname, char *usage_str) {
  using namespace std;
  string flags = "[";
  int i = 0;
  for(i=0;exs[i].name != NULL;i++) {
    string s = string("--") + string(exs[i].name);
    flags += (s + ((exs[i+1].name != NULL) ? string("|") : string("]")));
  }
  std::cout << "usage: " << sysname << " " << flags << " " << usage_str 
      << std::endl;
   

  for(i=0;exs[i].name != NULL;i++) {
     string s = string("--") + string(exs[i].name);
     std::cout << "\t" << setw(20) << left << s  << "\t: " << exs[i].msg 
         << std::endl;
   }
}

} // namespace hogwild

} // namespace hazy
