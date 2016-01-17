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

#ifndef HAZY_HOGWILD_FRONTEND_UTIL_H
#define HAZY_HOGWILD_FRONTEND_UTIL_H

#include <getopt.h>

namespace hazy {
namespace hogwild {

//! Extends the definition of option in getopt to add in a message.
struct extended_option
{
  const char *name;
  int has_arg,*flag,val;
  const char *msg;
};

//! convert extended options to regular options for use in getopt.
option* convert_extended_options(const extended_option* exs);


//! printout a usage string that lists all the options.
/*! extended options are all default options from above.
 * sysname is typically argv[0]
 * usage str contains the usage of the mandatory no flag options.
 */
void print_usage(const extended_option *exs, char* sysname, char *usage_str);

} // namespace hogwild

} // namespace hazy
#endif
