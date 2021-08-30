#ifndef RATATOSK_SHARED_PAIR_HPP
#define RATATOSK_SHARED_PAIR_HPP

#include "PairID.hpp"

class SharedPairID {

	public:

        class SharedPairID_const_iterator : public std::iterator<std::forward_iterator_tag, uint32_t> {

            friend class SharedPairID;

            public:

                SharedPairID_const_iterator();
                SharedPairID_const_iterator(const SharedPairID_const_iterator& o);

                SharedPairID_const_iterator& operator=(const SharedPairID_const_iterator& o);

                uint32_t operator*() const;

                SharedPairID_const_iterator operator++(int);
                SharedPairID_const_iterator& operator++();

                bool operator==(const SharedPairID_const_iterator& o) const;
                bool operator!=(const SharedPairID_const_iterator& o) const;

                void clear();

            private:

            	SharedPairID_const_iterator(const SharedPairID* spid, const bool begin);

            	PairID::const_iterator global_pid_it, local_pid_it;

            	const PairID* global_pid;
                const PairID* local_pid;

            	bool l_valid, g_valid;
        };

        typedef SharedPairID_const_iterator const_iterator;

		SharedPairID();

		SharedPairID(const SharedPairID& o);
		SharedPairID(SharedPairID&& o);

        SharedPairID(PairID&& pid);
        SharedPairID(const PairID& pid);

		~SharedPairID();

		void clear();

		SharedPairID& operator=(const SharedPairID& o);
		SharedPairID& operator=(SharedPairID&& o);

        //SharedPairID& operator=(PairID&& pid);
        //SharedPairID& operator=(const PairID& pid);

		bool operator==(const SharedPairID& o) const;

        SharedPairID operator|(const SharedPairID& rhs) const;
        SharedPairID& operator|=(const SharedPairID& rhs);

        inline SharedPairID operator+(const SharedPairID& rhs) const {

            return this->operator|(rhs);
        }

        inline SharedPairID& operator+=(const SharedPairID& rhs){

            return this->operator|=(rhs);
        }

        SharedPairID operator|(const PairID& rhs) const;
        SharedPairID& operator|=(const PairID& rhs);

        inline SharedPairID operator+(const PairID& rhs) const {

            return this->operator|(rhs);
        }

        inline SharedPairID& operator+=(const PairID& rhs){

            return this->operator|=(rhs);
        }

		size_t getSizeInBytes() const;

		void add(const size_t id);
		bool contains(const size_t pair_id) const;

		size_t maximum() const;
		size_t minimum() const;

		size_t size() const;

		void runOptimize();

        PairID toPairID() const;
        vector<uint32_t> toVector() const;

        pair<const PairID*, const PairID*> getPairIDs() const;

        const PairID& getLocalPairID() const;

        void setGlobalPairID(const PairID* g_pid);
        void setGlobalPairID(const PairID& g_pid);

        void setLocalPairID(const PairID& pid);

        bool write(ostream& stream_out) const;
        bool read(istream& stream_in);

        void forceRoaringInternal();

		SharedPairID::const_iterator begin() const;
		SharedPairID::const_iterator end() const;

        inline size_t cardinality() const {

            return size();
        }

        inline bool isEmpty() const {

            return (cardinality() == 0);
        }

	private:

        inline bool hasGlobalOwnership() const {

            return static_cast<bool>((global_pid & flagMask) == ptr_owner);
        }

        inline const PairID* getGlobalPairID() const {

            return reinterpret_cast<const PairID*>(global_pid & ptrMask);
        }

        static inline size_t l_approximate_log2(size_t v) {

            static const uint8_t tab64[64] = {
                63,  0, 58,  1, 59, 47, 53,  2,
                60, 39, 48, 27, 54, 33, 42,  3,
                61, 51, 37, 40, 49, 18, 28, 20,
                55, 30, 34, 11, 43, 14, 22,  4,
                62, 57, 46, 52, 38, 26, 32, 41,
                50, 36, 17, 19, 29, 10, 13, 21,
                56, 45, 25, 31, 35, 16,  9, 12,
                44, 24, 15,  8, 23,  7,  6,  5
            };

            v |= v >> 1;
            v |= v >> 2;
            v |= v >> 4;
            v |= v >> 8;
            v |= v >> 16;
            v |= v >> 32;

            return tab64[(static_cast<size_t>((v - (v >> 1))*0x07EDD5E59A4E28C2)) >> 58];
        }


        static const size_t shiftMaskBits; // 1 bit

        static const uintptr_t ptr_owner; // Flag 1
        static const uintptr_t ptr_user; // Flag 0

        static const uintptr_t flagMask;
		static const uintptr_t ptrMask;

		uintptr_t global_pid;
		PairID local_pid;
};

#endif