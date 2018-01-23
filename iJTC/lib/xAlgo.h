
#ifndef xAlgo_H
#define xAlgo_H

namespace xAlgo {

	int binarySearch(float key, float* arr, int i_max, int i_min ){
		if(key<arr[i_min]) return -1;
		if(key>arr[i_max]) return -1;
		//		assert(key< arr[i_max] ||!"overflow!");
		//		assert(key> arr[i_min] ||!"underflow!");
		int mid = floor(float(i_max +i_min)/2);
		//	cout<<mid<<endl;
		if(mid == i_min ) return mid;
		if( arr[mid]> key) return binarySearch(key, arr, mid, i_min);
		else if( arr[mid] < key) return binarySearch(key, arr, i_max, mid);
		else return mid;
	}

	int binarySearch(float key, std::vector<float> & arr, int i_max, int i_min ){
		if(key<arr.at(i_min)) return -1;
		if(key>arr.at(i_max)) return -1;
		int mid = floor(float(i_max +i_min)/2);
		//	cout<<mid<<endl;
		if(mid == i_min ) return mid;
		if( arr[mid]> key) return binarySearch(key, arr, mid, i_min);
		else if( arr[mid] < key) return binarySearch(key, arr, i_max, mid);
		else return mid;
	}

	template <typename T> long BinarySearch(long n, const T *array, T value)
	{
		// 	Binary search in an array of n values to locate value.
		//
		// Array is supposed  to be sorted prior to this call.
		// If match is found, function returns position of element.
		// If no match found, function gives nearest element smaller than value.

		const T* pind;
		pind = std::lower_bound(array, array + n, value);
		if ( (pind != array + n) && (*pind == value) )
			return (pind - array);
		else
			return (pind - array - 1);
	}
}

#endif
