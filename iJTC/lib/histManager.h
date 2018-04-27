
class histManager {
		public: 
				histManager(){};
				template <typename T> T* regHist (const char* name, const char* title, int nbin, float x, float y);
				template <typename T> T* get(const char* name);
				void sumw2();
				void write();
		public:
				std::vector<TH1*> th1;
				std::map<const char*, TH1*> hkey;
};

template <typename T> T* histManager::regHist (const char* name, const char* title, int nbin, float x, float y){
		T *h = new T(name, title, nbin, x, y);
		th1.push_back(dynamic_cast<TH1*>(h));
		hkey[name] = th1.back();
		return h;
}
template <typename T> T* histManager::get (const char* name){
		return dynamic_cast<T*>(hkey[name]);
}

void histManager::sumw2(){
		for(auto x=th1.begin(); x!= th1.end(); ++x){
				(*x)->Sumw2();
		}
}
void histManager::write(){
		for(auto x=th1.begin(); x!= th1.end(); ++x) (*x)->Write();
}
