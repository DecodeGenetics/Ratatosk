#ifndef RATATOSK_TINYBF_HPP
#define RATATOSK_TINYBF_HPP

#include <iostream>

#include "libpopcnt.h"

template<typename T>
class TinyBloomFilter {

	public:

		TinyBloomFilter() : table(nullptr), sz_h(0) {}

		TinyBloomFilter(const size_t nb_elem, const size_t bits_per_elem) : table(nullptr), sz_h(0) {

		    if ((nb_elem != 0) && (bits_per_elem != 0)){

		        uint64_t nb_h = static_cast<uint64_t>((bits_per_elem * log(2)));

		        nb_h += static_cast<uint64_t>(fpp(bits_per_elem, nb_h) >= fpp(bits_per_elem, nb_h+1));

		        if (nb_h > 255){

		        	cerr << "TinyBloomFilter::TinyBloomFilter(): This class does not support using > 255 hash functions" << endl;
		        	clear();

		        	exit(1);
		        }

		        nb_h &= 0x00000000000000ffULL;
		        sz_h = static_cast<uint64_t>(max(rndup(bits_per_elem * nb_elem), static_cast<size_t>(64)));

		        if (sz_h > 0x0000ffffffffffff){

		        	cerr << "TinyBloomFilter::TinyBloomFilter(): Required size of " << sz_h << " bits cannot be handled by this class." << endl;
		        	clear();

		        	exit(1);
		        }

		        table = new uint64_t[sz_h / 64]();

		        sz_h = (sz_h << 8) | nb_h;
		    }
		}

		TinyBloomFilter(const TinyBloomFilter& o, const bool reset_bits = false) : table(nullptr), sz_h(o.sz_h) {

			if (o.table != nullptr){

				table = new uint64_t[(sz_h >> 8)/ 64]();

				if (!reset_bits) memcpy(table, o.table, ((sz_h >> 8) / 64) * sizeof(uint64_t));
			}
		}

		TinyBloomFilter(TinyBloomFilter&& o) : table(o.table), sz_h(o.sz_h) {

			o.table = nullptr;

			o.clear();
		}

		~TinyBloomFilter() {

			clear();
		}

		TinyBloomFilter<T>& operator=(const TinyBloomFilter& o){

			if (&o != this){

				clear();

				sz_h = o.sz_h;

				if (o.table != nullptr){

					table = new uint64_t[(sz_h >> 8) / 64];

					memcpy(table, o.table, ((sz_h >> 8) / 64) * sizeof(uint64_t));
				}
			}

			return *this;
		}

		TinyBloomFilter<T>& operator=(TinyBloomFilter&& o){

			if (&o != this){

				clear();

				table = o.table;
				sz_h = o.sz_h;

				o.table = nullptr;

				o.clear();
			}

			return *this;
		}

		inline void clear() {

			if (table != nullptr) delete[] table;

			table = nullptr;
			sz_h = 0;
		}

		inline void reset(){

			if (table != nullptr) memset(table, 0, ((sz_h >> 8) / 64) * sizeof(uint64_t));
		}

		void insert(const T& elem) {

			const uint64_t nb_h = sz_h & 0x00000000000000ffULL;
			const uint64_t mask = (sz_h >> 8) - 1;
			const uint64_t hv_2 = XXH64(&elem, sizeof(T), r_2);

			uint64_t i = 0;
			uint64_t hv_1 = XXH64(&elem, sizeof(T), r_1);

			for (; i != nb_h; ++i){

				if ((table[(hv_1 & mask) >> 6] & (1ULL << (hv_1 & 0x3fULL))) == 0) break;
				hv_1 += hv_2;
			}

			for (; i != nb_h; ++i){

				table[(hv_1 & mask) >> 6] |= 1ULL << (hv_1 & 0x3FULL);
				hv_1 += hv_2;
			}
		}

		bool query(const T& elem) const {

			const uint64_t nb_h = sz_h & 0x00000000000000ffULL;
			const uint64_t mask = (sz_h >> 8) - 1;
			const uint64_t hv_2 = XXH64(&elem, sizeof(T), r_2);

			uint64_t i = 0;
			uint64_t hv_1 = XXH64(&elem, sizeof(T), r_1);

			for (; i != nb_h; ++i){

				if ((table[(hv_1 & mask) >> 6] & (1ULL << (hv_1 & 0x3fULL))) == 0) break;
				hv_1 += hv_2;
			}

			return (i == nb_h);
		}

		uint64_t or_cardinality_bits(const TinyBloomFilter& o) const {

			if ((o.sz_h != sz_h) || (table == nullptr) || (o.table == nullptr))  return 0;

			const uint64_t sz_table = (sz_h >>  8) / 64;

			uint64_t pop = 0;

			for (uint64_t i = 0; i != sz_table; ++i){

				const uint64_t word = table[i] | o.table[i];

				pop += popcnt(&word, sizeof(uint64_t));
			}

			return pop;
		}

		size_t and_cardinality_bits(const TinyBloomFilter& o) const {

			if ((o.sz_h != sz_h) || (table == nullptr) || (o.table == nullptr))  return 0;

			const uint64_t sz_table = (sz_h >> 8) / 64;

			uint64_t pop = 0;

			for (uint64_t i = 0; i != sz_table; ++i){

				const uint64_t word = table[i] & o.table[i];

				pop += popcnt(&word, sizeof(uint64_t));
			}

			return pop;
		}

		inline bool isEmpty() const {

			bool isEmpty = (table != nullptr);

			for (uint64_t i = 0; (i < sz_h) && isEmpty; ++i) isEmpty = (table[i] != 0);

			return isEmpty;
		}

		inline size_t cardinality_bits() const {

			if (table == nullptr) return 0;

			return popcnt(table, ((sz_h >> 8) / 64) * sizeof(uint64_t));
		}

		inline size_t getNumberHashFunctions() const {

			return (sz_h >> 8);
		}

		size_t getSizeInBytes() const {

			return sizeof(uint64_t*) + sizeof(uint64_t) + (sz_h >> 8);
		}

		bool write(ostream& stream_out) const {

		    if (stream_out.good()) stream_out.write(reinterpret_cast<const char*>(&sz_h), sizeof(uint64_t));
		    if (stream_out.good() && ((sz_h >> 8) != 0)) stream_out.write(reinterpret_cast<const char*>(table), (sz_h >> 8) * sizeof(uint64_t));

		    return stream_out.good();
		}

		bool read(istream& stream_in) {

		    if (stream_in.good()){

		        clear();

		        if (stream_in.good()) stream_in.read(reinterpret_cast<char*>(&sz_h), sizeof(uint64_t));

		        if (stream_in.good() && ((sz_h >> 8) != 0)) {

		        	table = new uint64_t[(sz_h >> 8)];

		        	stream_in.read(reinterpret_cast<char*>(table), (sz_h >> 8) * sizeof(uint64_t));
		        }
		    }

		    return stream_in.good();
		}

	private:

        inline double fpp(const uint64_t bits, const uint64_t h) const {

        	const double bits_d = static_cast<double>(bits);
        	const double h_d = static_cast<double>(h);

            return pow(1 - exp(-(h_d / bits_d)), h_d);
        }

        uint64_t* table;

        uint64_t sz_h;

        static const uint64_t r_1;
        static const uint64_t r_2;
};

template<typename T> const uint64_t TinyBloomFilter<T>::r_1 = 49157;
template<typename T> const uint64_t TinyBloomFilter<T>::r_2 = 1610612741;

#endif