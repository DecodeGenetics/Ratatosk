#ifndef RATATOSK_PAIR_HPP
#define RATATOSK_PAIR_HPP

#include <random>

#include "roaring.hh"
#include "TinyBitmap.hpp"

class PairID {

    //Ensure that PairID::setPtrBmp is always allocated with an 8 bytes alignment
    struct alignas(8) Bitmap { Roaring r; };

    public:

        class PairID_const_iterator : public std::iterator<std::forward_iterator_tag, uint32_t> {

            friend class PairID;

            public:

                PairID_const_iterator();
                PairID_const_iterator(const PairID_const_iterator& o);

                ~PairID_const_iterator();

                void clear();

                PairID_const_iterator& operator=(const PairID_const_iterator& o);

                inline uint32_t operator*() const {

                    return static_cast<uint32_t>(ck_id);
                }


                PairID_const_iterator operator++(int);
                PairID_const_iterator& operator++();

                bool operator==(const PairID_const_iterator& o) const;
                bool operator!=(const PairID_const_iterator& o) const;

            private:

                const PairID* cs;

                size_t flag;

                size_t it_setBits;
                size_t cs_sz;

                uint64_t ck_id;

                const Roaring empty_roar;

                TinyBitmap t_bmp;

                Roaring::const_iterator it_roar;
                TinyBitmap::const_iterator it_t_bmp;

                PairID_const_iterator(const PairID* cs_, const bool beg);

                inline bool isInvalid() const {

                    return ((ck_id == 0xffffffffffffffff) || (it_setBits == cs_sz));
                }
        };

        typedef PairID_const_iterator const_iterator;

        PairID();
        PairID(const PairID& o); // Copy constructor
        PairID(PairID&& o); // Move  constructor

        ~PairID();

        PairID& operator=(const PairID& o);
        PairID& operator=(PairID&& o);

        bool operator==(const PairID& o) const;

        inline bool operator!=(const PairID& o) const {

            return !operator==(o);
        }

        PairID operator|(const PairID& rhs) const;
        PairID& operator|=(const PairID& rhs);

        PairID operator&(const PairID& rhs) const;
        PairID& operator&=(const PairID& rhs);

        PairID operator-(const PairID& rhs) const;
        PairID& operator-=(const PairID& rhs);

        inline PairID operator+(const PairID& rhs) const{

            return this->operator|(rhs);
        }

        inline PairID& operator+=(const PairID& rhs){

            return this->operator|=(rhs);
        }

        void clear();

        void add(const size_t pair_id);
        void remove(const size_t pair_id);
        bool contains(const size_t pair_id) const;

        inline bool isEmpty() const { return (size() == 0); }

        size_t maximum() const;
        size_t minimum() const;

        size_t size() const;

        size_t and_cardinality(const PairID& rhs) const;
        size_t and_cardinality(const PairID& rhs, const uint64_t min_shared) const;

        inline size_t cardinality() const {

            return size();
        }

        bool write(ostream& stream_out) const;
        bool read(istream& stream_in);

        size_t getSizeInBytes() const;

        const_iterator begin() const;
        const_iterator begin(const size_t id) const;
        const_iterator end() const;

        void runOptimize();

        Roaring toRoaring() const;
        vector<uint32_t> toVector() const;

        static PairID fastunion(const size_t sz, const PairID* const* p_id) {

            PairID pid_out;

            {
                vector<const Roaring*> v_roar_ptr;

                for (size_t i = 0; i < sz; ++i) {

                    if (p_id[i]->isBitmap()) v_roar_ptr.push_back(&(p_id[i]->getConstPtrBitmap()->r));
                }

                const size_t card_roaring = v_roar_ptr.size();

                if (card_roaring > 0) {

                    Bitmap* bitm = new Bitmap;

                    bitm->r = Roaring::fastunion(card_roaring, &(v_roar_ptr.front()));
                    pid_out.setBits = (reinterpret_cast<uintptr_t>(bitm) & pointerMask) | ptrBitmap;
                }
            }

            for (size_t i = 0; i < sz; ++i) {

                if (!p_id[i]->isBitmap()) pid_out |= *(p_id[i]);
            }

            return pid_out;
        }

        PairID forceRoaringInternal() const;
        PairID subsample(const uint64_t nb) const;

    private:

        void addSortedVector(const vector<uint32_t>& v);
        void removeSortedVector(const vector<uint32_t>& v);

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

        inline void releaseMemory(){

            const uintptr_t flag = setBits & flagMask;

            if (flag == ptrBitmap) delete getPtrBitmap();
            else if (flag == localTinyBitmap){

                uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
                TinyBitmap t_bmp(&setPtrTinyBmp);

                t_bmp.clear();
            }

            setBits = localBitVector;
        }

        inline void shrinkSize(){

            const uintptr_t flag = setBits & flagMask;

            if (flag == ptrBitmap) getPtrBitmap()->r.shrinkToFit();
            else if (flag == localTinyBitmap){

                uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
                TinyBitmap t_bmp(&setPtrTinyBmp);

                t_bmp.shrinkSize();

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            }
        }

        inline bool isBitmap() const { return ((setBits & flagMask) == ptrBitmap); }
        inline bool isTinyBitmap() const { return ((setBits & flagMask) == localTinyBitmap); }

        inline Bitmap* getPtrBitmap() const {

            return reinterpret_cast<Bitmap*>(setBits & pointerMask);
        }

        inline const Bitmap* getConstPtrBitmap() const {

            return reinterpret_cast<const Bitmap*>(setBits & pointerMask);
        }

        inline uint16_t* getPtrTinyBitmap() const {

            return reinterpret_cast<uint16_t*>(setBits & pointerMask);
        }

        static const size_t maxBitVectorIDs; // 64 bits - 3 bits for the color set type = 61
        static const size_t shiftMaskBits; // 3 bits

        // asBits and asPointer represent:
        // Flag 0 - A TinyBitmap which can contain up to 65488 uint
        // Flag 1 - A bit vector of 62 bits storing presence/absence of up to 62 integers
        // Flag 2 - A single integer
        // Flag 3 - A pointer to a CRoaring compressed bitmap which can contain up to 2^32 uint

        static const uintptr_t localTinyBitmap; // Flag 0
        static const uintptr_t localBitVector; // Flag 1
        static const uintptr_t localSingleInt; // Flag 2
        static const uintptr_t ptrBitmap; // Flag 3

        static const uintptr_t flagMask; // 0x7 (= 2^shiftMaskBits - 1)
        static const uintptr_t pointerMask; // 0xfffffffffffffff8 (= 2^64 - 1 - flagMask)

        uintptr_t setBits;
};

#endif
