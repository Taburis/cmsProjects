
class histManager {
		public: 
				histManager(){};
				template <typename T> T* regHist (const char* name, const char* title, int nbin, double x, double y);
				template <typename T> T* regHist (const char* name, const char* title, int nbin, const float *bins);
				template <typename T> T* regHist (const char* name, const char* title, int nx, double x1, double x2, 
								int ny, double y1, double y2);
				template <typename T> T* get(const char* name);
				template <typename T> int fill(const char * name, double x, double w=1){
						return get<T>(name)->Fill(x, w);
				}
				void sumw2();
				void write();
		public:
				std::vector<TH1*> th1;
				std::map<const char*, TH1*> hkey;
};

template <typename T> T* histManager::regHist (const char* name, const char* title, int nbin, double x, double y){
		T *h = new T(name, title, nbin, x, y);
		th1.push_back(dynamic_cast<TH1*>(h));
		hkey[name] = th1.back();
		return h;
}

template <typename T> T* histManager::regHist (const char* name, const char* title, int nbin, const float* bins){
		T *h = new T(name, title, nbin, bins);
		th1.push_back(dynamic_cast<TH1*>(h));
		hkey[name] = th1.back();
		return h;
}

template <typename T> T* histManager::regHist (const char* name, const char* title, int nx, double x1, double x2,
				int ny, double y1, double y2){
		T *h = new T(name, title, nx, x1, x2, ny, y1, y2);
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
