
#ifndef xPSet_cpp
#define xPSet_cpp
#include "any.hpp"
#include <unordered_map>

class xPSet{
		public: 
				xPSet(){}
				bool exists(std::string name){
						return  table.find(name) == table.end() ? 0 : 1;
				}
				template <typename T>
						bool addPara(const char* pname, const T & val){
								if(exists(pname)){
										std::cout<<"parameter: "<<pname<<" exists in the Set"<<endl;
									   	return 0;
								}
								table[pname] = any(T(val));
								return 1;
						}
				template <typename T>
						bool setPara(const char* pname, const T & val){
								table[pname] = any(T(val));
								return 1;
						}
				template <typename T>
				T getPara(const char* pname){
						if(!exists(pname)) {
								std::cout<<"parameter '"<<pname<<"' not defined in the set!"<<endl;
								return 0;
						}
						return any_cast<T>(table[pname]);
				}

				template <typename T>
				T securePara(const char* pname, T def){
						if(!exists(pname)) {
								std::cout<<"parameter '"<<pname<<"' not defined! using defalut value."<<endl;
								setPara<T>(pname, def);
								return T(def);
						}
						return any_cast<T>(table[pname]);
				}
				xPSet* clone(){
						auto cps = new xPSet();
						cps->table = table;
						return cps;
				}
				xPSet& operator=(const xPSet & rhs){
						table = rhs.table;
						return (*this);
				};
		public: 
				std::string set_name;
				std::unordered_map<std::string, any> table;
};

#endif
