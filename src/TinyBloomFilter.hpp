#ifndef POLIFROST_TINYBF_HPP
#define POLIFROST_TINYBF_HPP

#include <iostream>

#include <bifrost/libpopcnt.h>

template<typename T>
class TinyBloomFilter {

	public:

		TinyBloomFilter() : table(nullptr), sz_table_bits(0), h(0), r_1(0), r_2(0), nb_elem_tot(0), nb_elem_curr(0) {}

		TinyBloomFilter(const size_t nb_elem, const size_t bits_per_elem) : table(nullptr), sz_table_bits(0), h(0), r_1(0), r_2(0), nb_elem_tot(0), nb_elem_curr(0) {

		    if ((nb_elem != 0) && (bits_per_elem != 0)){

			    /*std::random_device rd; //Seed
			    std::default_random_engine generator(rd()); //Random number generator
			    std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

			    r_1 = distribution(generator);
			    r_2 = distribution(generator);*/

			    r_1 = 49157;
			    r_2 = 1610612741;

		        sz_table_bits = max(rndup(bits_per_elem * nb_elem), static_cast<size_t>(64));

		        h = (int) (bits_per_elem * log(2));
		        h += (fpp(bits_per_elem, h) >= fpp(bits_per_elem, h+1));

		        table = new uint64_t[sz_table_bits / 64]();

			    nb_elem_tot = nb_elem;
		    }
		}

		TinyBloomFilter(const TinyBloomFilter& o, const bool reset_bits = false) : table(nullptr), sz_table_bits(o.sz_table_bits), h(o.h), r_1(o.r_1), r_2(o.r_2), nb_elem_tot(o.nb_elem_tot), nb_elem_curr(o.nb_elem_curr){

			if (o.table != nullptr){

				table = new uint64_t[sz_table_bits / 64]();

				if (!reset_bits) memcpy(table, o.table, (sz_table_bits / 64) * sizeof(uint64_t));
			}
		}

		TinyBloomFilter(TinyBloomFilter&& o) : table(o.table), sz_table_bits(o.sz_table_bits), h(o.h), r_1(o.r_1), r_2(o.r_2), nb_elem_tot(o.nb_elem_tot), nb_elem_curr(o.nb_elem_curr){

			o.table = nullptr;

			o.clear();
		}

		~TinyBloomFilter() {

			clear();
		}

		TinyBloomFilter<T>& operator=(const TinyBloomFilter& o){

			if (&o != this){

				clear();

				sz_table_bits = o.sz_table_bits;
				h = o.h;
				r_1 = o.r_1;
				r_2 = o.r_2;
				nb_elem_tot = o.nb_elem_tot;
				nb_elem_curr = o.nb_elem_curr;

				if (o.table != nullptr){

					table = new uint64_t[sz_table_bits / 64];

					memcpy(table, o.table, (sz_table_bits / 64) * sizeof(uint64_t));
				}
			}

			return *this;
		}

		TinyBloomFilter<T>& operator=(TinyBloomFilter&& o){

			if (&o != this){

				clear();

				sz_table_bits = o.sz_table_bits;
				h = o.h;
				r_1 = o.r_1;
				r_2 = o.r_2;
				table = o.table;
				nb_elem_tot = o.nb_elem_tot;
				nb_elem_curr = o.nb_elem_curr;

				o.table = nullptr;

				o.clear();
			}

			return *this;
		}

		inline void clear() {

			if (table != nullptr) delete[] table;

			table = nullptr;

			r_1 = 0;
			r_2 = 0;
			h = 0;
			sz_table_bits = 0;
			nb_elem_tot = 0;
			nb_elem_curr = 0;
		}

		inline void reset(){

			nb_elem_curr = 0;

			if (table != nullptr) memset(table, 0, (sz_table_bits / 64) * sizeof(uint64_t));
		}

		void insert(const T& elem) {

			int i = 0;

			uint64_t hv_1 = XXH64(&elem, sizeof(T), r_1);
			
			const uint64_t hv_2 = XXH64(&elem, sizeof(T), r_2);
			const uint64_t mask = sz_table_bits - 1;

			for (; i != h; ++i){

				if ((table[(hv_1 & mask) >> 6] & (1ULL << (hv_1 & 0x3fULL))) == 0) break;
				hv_1 += hv_2;
			}

			nb_elem_curr += (i != h);

			for (; i != h; ++i){

				table[(hv_1 & mask) >> 6] |= 1ULL << (hv_1 & 0x3FULL);
				hv_1 += hv_2;
			}
		}

		bool query(const T& elem) const {

			int i = 0;

			uint64_t hv_1 = XXH64(&elem, sizeof(T), r_1);
			
			const uint64_t hv_2 = XXH64(&elem, sizeof(T), r_2);
			const uint64_t mask = sz_table_bits - 1;

			for (; i != h; ++i){

				if ((table[(hv_1 & mask) >> 6] & (1ULL << (hv_1 & 0x3fULL))) == 0) break;
				hv_1 += hv_2;
			}

			return (i == h);
		}

		/*TinyBloomFilter union_bf(const TinyBloomFilter& o){

			if ((o.sz_table_bits != sz_table_bits) || (o.r_1 != r_1) || (o.r2 != r_2) || (o.h != h)) return TinyBloomFilter();

			const size_t sz_table = sz_table_bits / 64;

			for (size_t i = 0; i != sz_table; ++i) table[i] |= o.table[i];

			const size_t popcnt = popcount(table, sz_table * sizeof(uint64_t));

			nb_elem_curr = -sz_table_bits * log(1 - popcnt/sz_table_bits) / h;
		}*/

		size_t or_cardinality_bits(const TinyBloomFilter& o) const {

			if ((o.sz_table_bits != sz_table_bits) || (o.r_1 != r_1) || (o.r_2 != r_2) || (o.h != h)) return 0;
			if ((table == nullptr) || (o.table == nullptr))  return 0;

			const size_t sz_table = sz_table_bits / 64;

			size_t pop = 0;

			for (size_t i = 0; i != sz_table; ++i){

				const uint64_t word = table[i] | o.table[i];

				pop += popcnt(&word, sizeof(uint64_t));
			}

			return pop;
		}

		size_t and_cardinality_bits(const TinyBloomFilter& o) const {

			if ((o.sz_table_bits != sz_table_bits) || (o.r_1 != r_1) || (o.r_2 != r_2) || (o.h != h)) return 0;
			if ((table == nullptr) || (o.table == nullptr))  return 0;

			const size_t sz_table = sz_table_bits / 64;

			size_t pop = 0;

			for (size_t i = 0; i != sz_table; ++i){

				const uint64_t word = table[i] & o.table[i];

				pop += popcnt(&word, sizeof(uint64_t));
			}

			return pop;
		}

		inline bool isFull() const {

			return ((table != nullptr) && (nb_elem_curr >= nb_elem_tot));
		}

		inline bool isEmpty() const {

			return ((table == nullptr) || (nb_elem_curr == 0));
		}

		inline size_t cardinality() const {

			return nb_elem_curr;
		}

		inline size_t cardinality_bits() const {

			return popcnt(table, (sz_table_bits / 64) * sizeof(uint64_t));
		}

	private:

        inline double fpp(const size_t bits, const int h) const {

            return pow(1-exp(-((double)h)/((double)bits)),(double)h);
        }

        uint64_t* table;

        uint64_t sz_table_bits;
        uint64_t r_1;
        uint64_t r_2;

        size_t nb_elem_tot;
        size_t nb_elem_curr;

        int h;
};

#endif